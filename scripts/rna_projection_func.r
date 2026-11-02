library(dplyr)
library(tibble)
library(tidyr)

# TG packages
library(metacell)
library(anndata)
library(tgstat)
library(tgutil)
library(misha)

library(ggplot2)
library(ggrepel)
library(ggpubr)
library(cowplot)
library(zoo)
require(grid)

gsetroot(mm10_trackdb)

load_data <- function(seq = "mars", file_prefix) {
    metacells <- read_h5ad(sprintf("./data/MC_models/%s_%s_metacells_projected.h5ad", seq, file_prefix))
    cells <- read_h5ad(sprintf("./data/MC_models/%s_%s_cells.h5ad", seq, file_prefix))

    legc <- log2(t(data.matrix(metacells$X)) + 1e-5)
    lpgc <- log2(t(data.matrix(metacells$layers[["projected_fraction"]])) + 1e-5)
    title <- gsub("_", " ", file_prefix)
    mc_to_type <- setNames(metacells$obs$projected_type, rownames(metacells$obs))

    return(list(cells_X = cells$X, cells_obs = cells$obs, legc = legc, lpgc = lpgc, title = title, mc_to_type = mc_to_type))
}

# Load RNA_mod data, if reload_data is TRUE, the data will be reloaded from scratch
init_RNA_mod <- function(reload_data){
    
    if (reload_data){
        
        RNA_mod = list()
        RNA_mod$WT_4N <- load_data("mars", "WT_4N_ff")
        RNA_mod$DNMT_TKO_4N <- load_data("mars", "DNMT_TKO_4N_ff")
        RNA_mod$DNMT_DKO_4N <- load_data("mars", "DNMT_DKO_4N_ff")
        RNA_mod$DNMT_TKO_KO <- load_data("mars", "DNMT_TKO_KO_ff")
        RNA_mod$DNMT_TKO_host <- load_data("mars", "DNMT_TKO_host_ff")
        RNA_mod$DNMT_DKO_KO <- load_data("mars", "DNMT_DKO_KO_ff")
        RNA_mod$DNMT_DKO_host <- load_data("mars", "DNMT_DKO_host_ff")

        RNA_mod$multi.WT <- load_data("multiome", "WT")
        RNA_mod$multi.DNMT_TKO_KO <- load_data("multiome", "DNMT_TKO_KO")
        RNA_mod$multi.DNMT_DKO_KO <- load_data("multiome", "DNMT_DKO_KO")
        RNA_mod$multi.TET_TKO_KO <- load_data("multiome", "TET_TKO_KO")
        
        saveRDS(RNA_mod, file = out_path("./data/output/RNA_mod.RDS"))
        
    } else {
        RNA_mod = readRDS(file = in_path("./data/output/RNA_mod.RDS"))
    }

    return(RNA_mod)
}

load_gene_sets <- function(){

    scdb_init("./data/output/", force_reinit = T)
    gset_raw = scdb_gset("query_atlas")
    gset_raw = names(gset_raw@gene_set)

    TKO_kosher = gintervals.load("intervs.global.tss") %>% filter(!chrom %in% c("chr6", "chrX")) %>% pull(geneSymbol)
    DKO_kosher = gintervals.load("intervs.global.tss") %>% filter(!chrom %in% c("chr1", "chr14", "chrX")) %>% pull(geneSymbol)

    gset = gset_raw[gset_raw %in% intersect(TKO_kosher, DKO_kosher)]

    weber_genes = c("Sohlh2", "Tex13", "Rpl10l", "Dazl", "Asz1", "Hormad1", "Gtsf1", "Tex19.1", 
                "Topaz1", "Rnf17", "Piwil2", "Mov10l1", "Asz1", "Ddx4", "Mael", "Fkbp6", "Gpat2")

    imprinted_genes = c("Dnmt3a", "Dnmt3b", # control
                        "Nnat", "Peg3", "Phlda2", "Kcnq1ot1", "Dlk1", "Impact", "Snrpn",
                        "Airn", "Meg3", "Inpp5f", "Gnas", "Mest", "Nap1l5", "Cdkn1c", "Plagl1", "Peg10", "Ifg2r", "Rian")

    # Modules from TET paper
    epiblast_module = c("Dtx1", "Enpp3", "Foxd3", "L1td1", "Mal2", "Mkrn1", "Pdzd4", "Pim2", "Pipox", "Pou3f1", "Slc35d3", "Slc35f2", "Slc7a3", "Ugt8a", "Utf1")
    nm_module = c("Apln", "Cyp26a1", "Eomes", "Evx1", "Fgf8", "Frzb", "Fst", "Gas1", "Mesp1", "Mixl1", "Snai1", "Sp5", "T", "Tdgf1", "Wnt3")

    return(list("gset" = gset, "weber_genes" = weber_genes, "imprinted_genes" = imprinted_genes, 
                "epiblast_module" = epiblast_module, "nm_module" = nm_module))

}

# Downsampled (5000 UMIs) atlas matrix and atlas cell metadata (cell -> embryo, metacell, cell type).
# With re_downsample_mut_matrices = TRUE the per-dataset downsampled matrices are also rewritten to
# data/output/<dataset>_ds.RDS, which plot_density_plots reads.
init_ds_mats <- function(re_downsample_atlas_mat, re_downsample_mut_matrices){

    ref_atlas_path = paste0(spatwt_dir, "/")
    mc_wt = read_h5ad(paste0(ref_atlas_path, "spatwt.metacells.h5ad"))
    cells_wt = read_h5ad(paste0(ref_atlas_path, "spatwt.cells.h5ad"))

    if (re_downsample_atlas_mat){
        cells_wt_X = cells_wt$X
        mat_wt = t(data.matrix(cells_wt_X))

        a_umat = mat_wt
        a_ds3k = scm_downsamp(a_umat, 5e+3)
        a_ds3k = as.matrix(a_ds3k)

        saveRDS(a_ds3k, out_path("./data/output/atlas_ds.RDS"))
    } else {
        a_ds3k = readRDS(in_path("./data/output/atlas_ds.RDS"))
    }

    a_cl_to_mc = cells_wt$obs %>% select(embryo, metacell_name) %>% rename(metacell = metacell_name) %>% rownames_to_column("cell") %>%
            mutate(metacell = as.character(metacell))
    a_mc_to_ct = mc_wt$obs %>% select(projected_type) %>% 
            rename(cell_type = projected_type) %>% 
            rownames_to_column("metacell") %>% 
            mutate(cell_type = as.character(cell_type))

    a_md = full_join(a_cl_to_mc, a_mc_to_ct)
    a_md$cell_type[a_md$cell_type == "Early nascent mesoderm"] = "Nascent mesoderm Early"

    if (re_downsample_mut_matrices){
        for (name in names(RNA_mod)){
                subset_nm = name
                cells = RNA_mod[[name]]$cells_X
                umat = t(data.matrix(cells))

                # Generate downsampled single cell matrices
                ds3k = scm_downsamp(umat, 5e+3)
                ds3k = as.matrix(ds3k)

                saveRDS(ds3k, out_path(sprintf("./data/output/%s_ds.RDS", subset_nm)))
        }
    }

    return(list("a_ds3k" = a_ds3k, "a_md" = a_md))

}

create_density_plot <- function(name, ds3k, md, a_ds3k, a_md, module, query_color, x_limits, y_limits, this_cell_type, x_label, annotate_y_axis = T) {
    all_query_umis <- colSums(ds3k[module,])
    all_atlas_umis <- colSums(a_ds3k[module,])
    
    geno_cls <- md %>%
        filter(cell_type == this_cell_type & cell %in% names(all_query_umis)) %>%
        pull(cell)

    a_cls <- a_md %>%
        filter(cell_type == this_cell_type & cell %in% names(all_atlas_umis)) %>%
        pull(cell)

    ct_query_umis <- log2(all_query_umis[geno_cls] / 5000 + 2e-4)
    ct_atlas_umis <- log2(all_atlas_umis[a_cls] / 5000 + 2e-4)

    plot_df <- data.frame(
        origin = c(rep("query", length(ct_query_umis)), rep("atlas", length(ct_atlas_umis))),
        score = c(ct_query_umis, ct_atlas_umis)
    )
    
    p <- ggplot(plot_df, aes(x = score, linetype = origin)) +
        geom_density(bw = 0.5) +
        xlim(x_limits) + ylim(y_limits) +
        labs(x = NULL, y = NULL, title = gsub("_", " ", name)) +
        scale_linetype_manual(values = c("atlas" = 2, "query" = 1)) +
        theme_classic(base_size = 15)

    if (annotate_y_axis){
        p = p + theme(axis.text.y = element_text(size = 10), 
                      plot.title = element_text(size = 16, hjust = 0.5), legend.position = "none")
    }
    
    if (!annotate_y_axis){
        p = p + theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(), 
                      plot.title = element_text(size = 18, hjust = 0.5), legend.position = "none")
    }
    
    return(p)
}

plot_density_plots <- function(RNA_mod, gn_sets, ds_mats, annotate_y_axis = T){

    # The published panels score the modules restricted to the feature gene set (gn_sets$gset), which excludes
    # genes on the chromosomes that are trisomic in the TKO and DKO lines (chr6, chrX; chr1, chr14, chrX).
    # Without the filter the TKO 4N distribution shifts right, driven by Pim2, Pdzd4, Slc7a3 (chrX) and Mkrn1 (chr6).
    epiblast_module = gn_sets$epiblast_module[gn_sets$epiblast_module %in% gn_sets$gset]
    nm_module = gn_sets$nm_module[gn_sets$nm_module %in% gn_sets$gset]

    epi_plots = list()
    nm_plots = list()
    
    for (name in c('DNMT_TKO_4N', 'DNMT_DKO_4N', 'DNMT_TKO_KO', 'DNMT_TKO_host')){

        cells = RNA_mod[[name]][["cells_obs"]]
        mcells = RNA_mod[[name]][["mcells"]]
        color = RNA_mod[[name]][["color"]]
        ds3k = readRDS(in_path(sprintf("./data/output/%s_ds.RDS", name)))

        # Associate each cell with a cell type
        cl_to_mc = cells %>% select(metacell_name, genotype) %>% rename(metacell = metacell_name) %>% rownames_to_column("cell") %>%
                    mutate(metacell = as.character(metacell))
        
        
        mc_to_ct = RNA_mod[[name]]$mc_to_type %>% enframe(name = "metacell", value = "cell_type")
        
        md = full_join(cl_to_mc, mc_to_ct)

        p_epi <- create_density_plot(name, ds3k, md, ds_mats$a_ds3k, ds_mats$a_md, epiblast_module, color, c(-10, -6), c(0, 1), "Epiblast", "Epiblast", annotate_y_axis = annotate_y_axis)
        p_nm <- create_density_plot(name, ds3k, md, ds_mats$a_ds3k, ds_mats$a_md, nm_module, color, c(-10, -6), c(0, 1), "Nascent mesoderm Early", "ENM", annotate_y_axis = annotate_y_axis)
        
        epi_plots[[name]] = p_epi
        nm_plots[[name]] = p_nm
    }

    figure_epi = ggarrange(plotlist = epi_plots, nrow = 1)
    figure_epi = annotate_figure(figure_epi, left = textGrob("Density", rot = 90, vjust = 1, gp = gpar(cex = 1.5)),
                        bottom = textGrob("Score", gp = gpar(cex = 1.5)))

    figure_nm = ggarrange(plotlist = nm_plots, nrow = 1)
    figure_nm = annotate_figure(figure_nm, left = textGrob("Density", rot = 90, vjust = 1, gp = gpar(cex = 1.5)),
                        bottom = textGrob("Score", gp = gpar(cex = 1.5)))

    show_fig("./figs/Fig2/epi_density.png", figure_epi, width = 10, height = 2.5)
    show_fig("./figs/Fig2/nm_density.png", figure_nm, width = 10, height = 2.5)
    invisible(list(epi = figure_epi, nm = figure_nm))
}

find_best_matches <- function(name, this_cell_type, gset){
    gc_per_embryo = read_h5ad(in_path(sprintf("./data/gc_per/gc_per_embryo_%s.h5ad", name)))
    gc_per_wt_embryo = read_h5ad(in_path("./data/gc_per/gc_per_embryo_atlas.h5ad"))

    # Filter gset
    gset = gset[gset %in% colnames(gc_per_embryo$X)]

    # Identify cell type columns and generate matrices
    obs_subset = gc_per_embryo$obs %>% filter(cell_type == this_cell_type) %>% rownames_to_column("index")
    wt_obs_subset = gc_per_wt_embryo$obs %>% filter(cell_type == this_cell_type & !embryo %in% c("", "NA", "empty")) %>% rownames_to_column("index")

    gc_per_embryo_mat = t(data.matrix(gc_per_embryo$X))
    gc_per_embryo_mat_subset = gc_per_embryo_mat[,obs_subset$index]
    colnames(gc_per_embryo_mat_subset) = obs_subset$embryo

    gc_per_wt_embryo_mat = t(data.matrix(gc_per_wt_embryo$X))
    gc_per_wt_embryo_mat_subset = gc_per_wt_embryo_mat[,wt_obs_subset$index]
    colnames(gc_per_wt_embryo_mat_subset) = wt_obs_subset$embryo

    legc_query = log2(gc_per_embryo_mat_subset + 1e-5)
    legc_wt = log2(gc_per_wt_embryo_mat_subset + 1e-5)

    # Calculate query~atlas correlation
    query_corrs = tgs_cor(legc_query[gset,], legc_wt[gset,])
    ordered_rows <- order(apply(query_corrs, 1, max), decreasing = TRUE)
    query_corrs <- query_corrs[ordered_rows,]
    
    best_matches <- c()
    best_query_corrs <- c()
    # Iterate over each row of the correlation matrix to choose as many best matches as query embryos
    for (i in 1:nrow(query_corrs)) {
        max_col <- names(which.max(query_corrs[i,]))
        best_matches <- c(best_matches, max_col)
        best_query_corrs <- c(best_query_corrs, max(query_corrs[i,]))
        query_corrs <- query_corrs[,colnames(query_corrs) != max_col]
    }
    
    # Calculate atlas~atlas correlation
    wt_corrs = tgs_cor(legc_wt[gset, best_matches], legc_wt[gset,])
    for (emb in best_matches){
        wt_corrs[emb, emb] = 0
    }
    ordered_rows <- order(apply(wt_corrs, 1, max), decreasing = TRUE)
    wt_corrs <- wt_corrs[ordered_rows,]
    
    best_wt_corrs <- c()
    # Iterate over each row of the correlation matrix to choose as many best matches as query embryos
    for (i in 1:nrow(wt_corrs)) {
        max_col <- names(which.max(wt_corrs[i,]))
        best_wt_corrs <- c(best_wt_corrs, max(wt_corrs[i,]))
        wt_corrs <- wt_corrs[,colnames(wt_corrs) != max_col]
    }
    
    return(list("query_embryos" = colnames(legc_query), "best_matches" = best_matches, 
                "best_query_corrs" = best_query_corrs, "best_wt_corrs" = best_wt_corrs))
    
}

plot_best_corrs <- function(name, best_matches_l){
    
    best_matches = best_matches_l[[name]]
    error_bars = data.frame("group" = c(rep("query", length(best_matches$query_embryos)), rep("control", length(best_matches$query_embryos))),
                            "corr" = c(1-best_matches$best_query_corrs, 1-best_matches$best_wt_corrs)) %>%
        group_by(group) %>% summarize(median = median(corr, na.rm = T),
                                  sd_pos = median + sd(corr, na.rm = T),
                                  sd_neg = median - sd(corr, na.rm = T))

    plot_df <- rbind(data.frame("corr" = 1 - best_matches$best_query_corrs, "group" = "query"), 
          data.frame("corr" = 1 - best_matches$best_wt_corrs, "group" = "control")) %>%
        mutate(group = factor(group, c("query", "control")))
    
    
    pval <- wilcox.test(best_matches$best_query_corrs, best_matches$best_wt_corrs, alternative = "two.sided")$p.value
    
    p <- ggplot() +
        geom_jitter(data = plot_df, aes(x = group, y = corr), position = position_jitter(width = 0.2, seed = 1)) +  # seed 1; the published panels used unseeded jitter
        geom_errorbar(data = error_bars, aes(x = group, ymin=sd_neg, ymax=sd_pos), width = 0.3) + 
        geom_errorbar(data = error_bars, aes(x = group, ymin = median, ymax = median), color = "black", size = 1.5, width = 0.7) + 
        labs(title = gsub("_", " ", name), subtitle = paste("pv =", formatC(pval, format = "e", digits = 1)), x = "", y = "") + ylim(0, 0.1) +
        theme_classic(base_size = 15) + theme(axis.text.y = element_text(size = 10), plot.title = element_text(size = 16, hjust = 0.5)) +
        scale_fill_identity()
    
    return(p)
}

plot_divergence_plots <- function(gset){
    best_matches_epi = lapply(c("DNMT_TKO_4N", "DNMT_DKO_4N", "DNMT_TKO_KO", "DNMT_TKO_host"), function(x) {
    find_best_matches(name = x, this_cell_type = "Epiblast", gset = gset)})
    names(best_matches_epi) = c("DNMT_TKO_4N", "DNMT_DKO_4N", "DNMT_TKO_KO", "DNMT_TKO_host")

    best_matches_nm = lapply(c("DNMT_TKO_4N", "DNMT_DKO_4N", "DNMT_TKO_KO", "DNMT_TKO_host"), function(x) {
        find_best_matches(name = x, this_cell_type = "Nascent mesoderm Early", gset = gset)})
    names(best_matches_nm) = c("DNMT_TKO_4N", "DNMT_DKO_4N", "DNMT_TKO_KO", "DNMT_TKO_host")

    epi_plots = list()
    nm_plots = list()

    for (name in names(best_matches_epi)){
        p <- plot_best_corrs(name, best_matches_epi)
        epi_plots[[name]] = p
    }

    for (name in names(best_matches_nm)){
        p <- plot_best_corrs(name, best_matches_nm)
        nm_plots[[name]] = p
    }

    options(repr.plot.width = 10, repr.plot.height = 2.5)
    figure_epi = ggarrange(plotlist = epi_plots, nrow = 1)
    figure_epi = annotate_figure(figure_epi, left = textGrob("Divergence score", rot = 90, vjust = 1, gp = gpar(cex = 1.5)))

    figure_nm = ggarrange(plotlist = nm_plots, nrow = 1)
    figure_nm = annotate_figure(figure_nm, left = textGrob("Divergence score", rot = 90, vjust = 1, gp = gpar(cex = 1.5)))

    show_fig("./figs/Fig2/epi_divergence.png", figure_epi, width = 10, height = 2.5)
    show_fig("./figs/Fig2/nm_divergence.png", figure_nm, width = 10, height = 2.5)
    invisible(list(epi = figure_epi, nm = figure_nm))
}

gen_obs_vs_proj_df <- function(gn_list, legc, lpgc, title, color){

    gn_list = gn_list[gn_list %in% rownames(legc)]

    plot_df_l = list()
    for (gn in gn_list){
        plot_df_s = data.frame("Gene" = gn, "obs" = unname(legc[gn,]), "proj" = unname(lpgc[gn,]))
        plot_df_l[[gn]] = plot_df_s
    }
            
    plot_df = do.call(rbind, plot_df_l) %>%
        pivot_longer(-Gene, names_to = "group", values_to = "expression") %>%
        group_by(Gene, group) %>% summarize(median_expression = median(expression)) %>% 
        mutate(facet = title, fill = ifelse(group == "obs", "deepskyblue", "grey20"), fill = as.character(fill))
                
    return(plot_df)
    
}

plot_obs_vs_proj_barplot <- function(gn_list, facets){
    
    names = gsub(" ", "_", facets)
    
    plot_df_l = lapply(names, function(x) {
        entry = RNA_mod[[x]]
        plot_df = gen_obs_vs_proj_df(gn_list, entry$legc, entry$lpgc, gsub("_", " ", x), entry$color)
        return(plot_df)
        }
    )
    
    full_plot_df = do.call(rbind, plot_df_l)
    full_plot_df$facet = factor(full_plot_df$facet, facets)
        
    gn_order_df = full_plot_df %>% filter(facet == "DNMT TKO 4N", group == "obs") %>% 
        arrange(median_expression)
    full_plot_df$Gene = factor(full_plot_df$Gene, gn_order_df$Gene) 
    full_plot_df = full_plot_df %>% mutate(median_expression = median_expression - log2(1e-5))
    
    p <- ggplot(full_plot_df, aes(x = median_expression, 
                         y = Gene, 
                         fill = fill)) + 
        geom_bar(stat = "identity", position = "dodge", width = 0.75) + 
        xlab("median(log2(expression))") + ylab("") + theme_minimal(base_size = 20) +
        scale_x_continuous(limits = c(0, -14 - log2(1e-5)), labels = function(x) {round(x + log2(1e-5), digits = 1)}) +
        facet_wrap(~facet, nrow = 1) + scale_fill_identity() +
        theme(legend.position = "None", axis.text.x = element_text(size = 14), plot.title = element_text(face = "bold", hjust = 0.5),
              plot.title.position = "plot", plot.margin = margin(t = 20, r = 10, b = 10, l = 10))
    
    return(p)
}

gen_diff_df <- function(gn_list, legc, lpgc, title, color){

    gn_list = gn_list[gn_list %in% rownames(legc)]

    plot_df_l = list()
    for (gn in gn_list){
        plot_df_s = data.frame("Gene" = gn, "obs" = unname(legc[gn,]), "proj" = unname(lpgc[gn,]))
        plot_df_l[[gn]] = plot_df_s
    }
    
    plot_df_l[[""]] = data.frame("Gene" = "", "obs" = NA, "proj" = NA)
            
    plot_df = do.call(rbind, plot_df_l) %>%
        mutate(diff = obs - proj) %>% select(-c(obs, proj)) %>%
        group_by(Gene) %>% summarize(median_diff = median(diff)) %>% 
        mutate(facet = title)
        
    return(plot_df)
    
}

plot_diff_barplot <- function(gn_list, facets){
    
    names = gsub(" ", "_", facets)
    
    plot_df_l = lapply(names, function(x) {
        entry = RNA_mod[[x]]
        plot_df = gen_diff_df(gn_list, entry$legc, entry$lpgc, gsub("_", " ", x), entry$color)
        return(plot_df)
        }
    )

    full_plot_df = do.call(rbind, plot_df_l)
    full_plot_df$facet = factor(full_plot_df$facet, facets)
    
    gn_order_df = full_plot_df %>% filter(facet == "DNMT TKO 4N", !Gene %in% c("", "Dnmt3a", "Dnmt3b")) %>% arrange(median_diff)
    full_plot_df = full_plot_df %>%
        mutate(Gene = factor(Gene, c(gn_order_df$Gene, "", "Dnmt3b", "Dnmt3a")),
               fill = "deepskyblue3")
    
    p <- ggplot(full_plot_df, aes(x = median_diff, 
                                  y = Gene, 
                                  fill = fill)) + 
        geom_bar(stat = "identity", position = "dodge", width = 0.75) + 
        xlab("median(log2(obs/proj))") + ylab("") + theme_minimal(base_size = 20) +
        xlim(-5,3) +
        facet_wrap(~facet, nrow = 1) + scale_fill_identity() +
        theme(legend.position = "None", axis.text.x = element_text(size = 14), plot.title = element_text(face = "bold", hjust = 0.5),
              plot.title.position = "plot", plot.margin = margin(t = 20, r = 10, b = 10, l = 10))
    
    return(p)
}
