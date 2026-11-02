library(dplyr)
library(tibble)
library(tidyr)

# TG packages
library(metacell)
library(anndata)
library(tgstat)
library(tgutil)

library(ggplot2)
library(ggpubr)
library(cowplot)
library(zoo)

library(data.table)
library(cli)
library(colorspace)

# Load cell type md
load_cell_type_md <- function(){
    color_key = fread(paste0(spatwt_dir, "/cell_type_colors.tsv"))
    color_key = rbind(color_key, data.frame("cell_type" = c("Early nascent mesoderm", "Late nascent mesoderm"), "color" = c("#C594BF", "#DFCDE4")))

    cell_type_ref = read.csv("./data/ref/cell_type_ref.csv")
    lineage_order = rev(c("ICM", "PGC", "Ectoderm", "Epiblast", "PS", "Nacent mesoderm", "Embryonic mesoderm", "ExE mesoderm", 
                        "Haematoendothelial progenitors", "Endothelial cells", "Blood", "Endoderm"))

    lineage_color_order = cell_type_ref %>% select(lineage, lineage_color) %>% distinct() %>% filter(lineage %in% lineage_order) %>% deframe()
    lineage_color_order = unname(lineage_color_order[lineage_order])

    cell_type_ref$lineage_color = factor(cell_type_ref$lineage_color, levels = lineage_color_order)
    cell_type_ref = cell_type_ref %>% filter(lineage != "Blood")

    emb_classification = fread("./data/ref/embryo_classification.csv", header = T) %>% select(1,2)

    return(list(color_key = color_key, cell_type_ref = cell_type_ref, emb_classification = emb_classification))
}

# Functions for plotting cell type composition of the models.

# Save data to data/composition, pdf and png files to figs/Fig1.
# Return embryo order for the plots.
plot_composition <- function(name, title, embryo_orders, ct_md = cell_type_md, annotate_dissection_time = TRUE, min_N = 100, width = 10, height = 5){

    # Import data
    if (name == "atlas"){
        ref_atlas_path = paste0(spatwt_dir, "/")
        mcells = read_h5ad(paste0(ref_atlas_path, "spatwt.metacells.h5ad"))
        cells = read_h5ad(paste0(ref_atlas_path, "spatwt.cells.h5ad"))
        
    } else {
        mcells = read_h5ad(sprintf("./data/MC_models/mars_%s_metacells_projected.h5ad", name))
        cells = read_h5ad(sprintf("./data/MC_models/mars_%s_cells.h5ad", name))
    }
    
    cells_obs = cells$obs
    
    # Prepare data for plotting
    plot_df = cells_obs %>% select(metacell_name, embryo) %>% rownames_to_column("id") %>%
        left_join(mcells$obs %>% rownames_to_column("metacell_name") %>% select(metacell_name, projected_type)) %>%
        rename(cell_type = projected_type) %>% left_join(ct_md$color_key) %>% group_by(embryo) %>% count(color) %>% left_join(ct_md$color_key)
    
    pre_f_cts = unique(plot_df$cell_type)
    plot_df = plot_df %>%
        filter(color %in% ct_md$cell_type_ref$color) %>%
        mutate(color = factor(color, levels = unique(ct_md$cell_type_ref$color)))
        
    # Filter out embryos with less than min_N cells
    pre_f_num = unique(plot_df$embryo) %>% length()
    enough_cls_embryos = plot_df %>% group_by(embryo) %>% summarize(N = sum(n)) %>% filter(N >= min_N) %>% pull(embryo) 
    plot_df = plot_df %>% filter(embryo %in% enough_cls_embryos)
    post_f_num = unique(plot_df$embryo) %>% length()
    cli_alert_info(text = sprintf("Filtered out %s embryos with less than %s cells", pre_f_num - post_f_num, min_N))
    
    plot_df = plot_df[,c("embryo", "cell_type", "color", "n")]
    saveRDS(plot_df, out_path(sprintf("./data/output/composition/%s.RDS", gsub("_ff", "", name))))
    
    # Manually set embryo order for DNMT_DKO_KO
    if (name == "DNMT_DKO_KO"){
        message("Setting manual order.")
        embryo_order = paste0("2NDKO_e", c(3, 2, 1, 6, 5, 10, 8, 7, 4))

    # Order embryo by the proportion of ExE mesoderm to early cell types
    } else {
        embryo_order = plot_df %>% 
            group_by(embryo) %>% mutate(f = n / sum(n)) %>%
            filter(cell_type %in% c("ExE mesoderm", "Amnion/Chorion progenitor", "Amnion/Chorion", "Allantois", 
                                    "Epiblast", "Primitive Streak", "Rostral neural plate", 
                                    "Early nascen mesoderm", "Nascent mesoderm Early", "Late nascent mesoderm", "Nascent mesoderm late")) %>%
            mutate(ct_group = ifelse(cell_type %in% c("ExE mesoderm", "Amnion/Chorion progenitor", "Amnion/Chorion", "Allantois"), "exe_meso", "early")) %>% 
            group_by(embryo, ct_group) %>%
            summarize(sum = sum(f)) %>% ungroup() %>%
            pivot_wider(names_from = ct_group, values_from = sum) %>% 
            replace(is.na(.), 0) %>%
            mutate(exe_meso_f = exe_meso + 0.01 / early + 0.01) %>% arrange(exe_meso_f) %>% pull(embryo) %>% as.character()
    }
    
    embryo_orders[[name]] = embryo_order
        
    plot_df = plot_df %>% 
        left_join(ct_md$cell_type_ref) %>% left_join(ct_md$emb_classification) %>%
        mutate(embryo = factor(embryo, levels = embryo_order)) %>% filter(!is.na(embryo))
    
    pdf_name = name

    annotation_colors = c("E7.5" = "grey70", "E8.5" = "grey30")
    annotation_plot <- ggplot(plot_df, aes(x = embryo, y = 1, fill = dissection_day)) +
        geom_tile() +
        scale_fill_manual(values = annotation_colors) +
        theme_void() + 
        theme(legend.position = "none")
        
    bar_plot <- ggplot(arrange(plot_df, lineage_color), aes(embryo, n)) + geom_bar(position = "fill", stat="identity", aes(fill=lineage_color)) +
        ggtitle(title) +  ylab("") +
        scale_fill_identity() + theme_void(base_size = 20) + theme(# axis.text.x=element_text(size = 10, angle = 90, vjust = 0.5, hjust=1), 
                                                                   axis.text.x=element_blank(),
                                                                   axis.ticks.x=element_blank(),
                                                                   plot.title = element_text(size = 20, hjust = 0.5, face = "bold"))
    
    if (annotate_dissection_time){
        library(patchwork)
        p <- bar_plot / annotation_plot + plot_layout(heights = c(1, 0.1))
    } else {
        p <- bar_plot
    }
    
    show_fig(sprintf("./figs/Fig1/%s.pdf", pdf_name), p, width = width, height = height)
    
    return(embryo_order)
}

# Plot all MARS composition plots, return embryo orders for the plots.
plot_mars_compositions <- function(cell_type_md){

    embryo_orders = list()
    # Single 2N KOs
    embryo_orders[["DNMT1"]] <- plot_composition("DNMT1_KO_ff", "DNMT1", embryo_orders = embryo_orders)
    embryo_orders[["UHRF1"]] <- plot_composition("UHRF1_KO_ff", "UHRF1", embryo_orders = embryo_orders)
    embryo_orders[["DNMT3A"]] <- plot_composition("DNMT3A_KO_ff", "DNMT3A", embryo_orders = embryo_orders)
    embryo_orders[["DNMT3B"]] <- plot_composition("DNMT3B_KO_ff", "DNMT3B", embryo_orders = embryo_orders)

    # 2N KOs
    embryo_orders[["DNMT_TKO"]] <- plot_composition("DNMT_TKO_KO_ff", "DNMT TKO", embryo_orders = embryo_orders)
    embryo_orders[["DNMT_DKO"]] <- plot_composition("DNMT_DKO_KO_ff", "DNMT DKO", embryo_orders = embryo_orders)

    # 4N KOs
    embryo_orders[["WT_4N"]] <- plot_composition("WT_4N_ff", "4N control", embryo_orders = embryo_orders)
    embryo_orders[["DNMT_TKO_4N"]] <- plot_composition("DNMT_TKO_4N_ff", "4N DNMT TKO", embryo_orders = embryo_orders)
    embryo_orders[["DNMT_DKO_4N"]] <- plot_composition("DNMT_DKO_4N_ff", "4N DNMT DKO", embryo_orders = embryo_orders)

    embryo_orders[["atlas"]] <- plot_composition("atlas", "atlas", embryo_orders = embryo_orders)

    return(embryo_orders)
}

organize_lineage_stats <- function(lin, ct_md = cell_type_md){
    
    facet_cts = ct_md$cell_type_ref %>% filter(lineage == lin) %>% pull(cell_type)
    facet_col = ct_md$cell_type_ref %>% filter(lineage == lin) %>% pull(lineage_color) %>% unique() %>% as.character()

    facet_col <- lighten(facet_col)
    facet_col2 <- darken(facet_col, amount = 0.2)

    genotypes = c("WT_4N", "DNMT_TKO_4N", "DNMT_DKO_4N", 
                  "DNMT_TKO_host", "DNMT_TKO_KO", "DNMT_DKO_host", "DNMT_DKO_KO",
                  "DNMT1_KO", "DNMT1_host", "UHRF1_KO", "UHRF1_host", "DNMT3A_KO", "DNMT3A_host", "DNMT3B_KO", "DNMT3B_host")
    titles = c("4N control", "4N TKO", "4N DKO", 
               "TKO host", "2N TKO", "DKO host", "2N DKO",
               "DNMT1 KO", "DNMT1 host", "UHRF1 KO", "UHRF1 host", "DNMT3A KO", "DNMT3A host", "DNMT3B KO", "DNMT3B host")
    names(titles) = genotypes

    composition_plot = list()
        for (genotype in genotypes){
            composition = readRDS(sprintf("./data/output/composition/%s.RDS", genotype))

            cts = facet_cts

            composition_mod = composition %>% filter(color %in% ct_md$cell_type_ref$color) %>%
                                  group_by(embryo) %>% select(-color) %>%
                                  mutate(f = n / sum(n)) %>%
                                  ungroup() %>% filter(cell_type %in% cts) 

        for (this_embryo in unique(composition$embryo)){
            existing_cts <- composition_mod %>% select(embryo, cell_type, f) %>% 
                                                      filter(f != 0) %>%
                                                      filter(embryo == this_embryo)

            missing_cts = setdiff(cts, existing_cts)
            if (length(missing_cts) >= 1){
                    missing_df <- data.frame("embryo" = this_embryo, 
                                             "cell_type" = missing_cts, 
                                             "n" = 0, 
                                             "f" = 0)
                    colnames(missing_df) = c("embryo", "cell_type", "n", "f")
                    composition_mod = rbind(composition_mod[,colnames(missing_df)], missing_df)                
            }
        }

        composition_mod = composition_mod %>%
                                group_by(embryo) %>% mutate(total_f = sum(f)) %>%
                                ungroup() %>%
                                select(-c(cell_type, n, f)) %>% distinct()

        embryos = embryo_orders[[gsub("_KO|_host", "", genotype)]]
        emb_to_diss = ct_md$emb_classification %>% filter(embryo %in% embryos) %>% deframe()

        composition_mod = composition_mod %>% left_join(enframe(emb_to_diss, name = "embryo", "diss_time")) %>% mutate(genotype = genotype)
        composition_plot[[genotype]] = composition_mod
    }

    plot_df = do.call(rbind, composition_plot) %>% mutate(title = titles[genotype], title = factor(title, titles))
    error_bars = plot_df %>% group_by(title, diss_time) %>% summarize(median = median(total_f, na.rm = T),
                                                                     sd_pos = median + sd(total_f, na.rm = T),
                                                                     sd_neg = median - sd(total_f, na.rm = T))
    
    return(list("plot_df" = plot_df, "error_bars" = error_bars, "facet_col" = facet_col, "facet_col2" = facet_col2))
}

produce_4N_violins <- function(lin, embryo_orders, ct_md = cell_type_md){
    
    stats = organize_lineage_stats(lin)
    
    plot_df = stats$plot_df
    error_bars = stats$error_bars
    
    plot_df_4N = plot_df %>% filter(title %in% c("4N control", "4N TKO", "4N DKO"))
    error_bars_4N = error_bars %>% filter(title %in% c("4N control", "4N TKO", "4N DKO"))
    
    pvals_4N <- plot_df_4N %>%
        group_by(title, diss_time) %>%
        summarize(
            p_value = wilcox.test(total_f, plot_df$total_f[plot_df$title == "4N control" & plot_df$diss_time == diss_time])$p.value
        )

        # Step 2: Add significance levels based on p-values
        pvals_4N <- pvals_4N %>%
            mutate(significance = case_when(
                p_value < 0.001 ~ "***",
                p_value < 0.01 ~ "**",
                p_value < 0.05 ~ "*",
                TRUE ~ "ns"
                )) %>%
        filter(title != "4N control")
    
    for (time in c("E7.5", "E8.5")){
        ylim = (plot_df_4N %>% filter(diss_time == time) %>% pull(total_f) %>% max) + 0.1
        p <- ggplot() +
                geom_errorbar(data = subset(error_bars_4N, diss_time == time), 
                              aes(x = title, ymin=sd_neg, ymax=sd_pos, color = diss_time), 
                              width = 0.3, position = position_dodge(width = 0.8)) + 
                geom_errorbar(data = subset(error_bars_4N, diss_time == time), 
                              aes(x = title, ymin = median, ymax = median, color = diss_time), 
                              size = 1.5, width = 0.5, position = position_dodge(width = 0.8)) + 
                geom_violin(data = subset(plot_df_4N, diss_time == time), 
                           aes(title, total_f, fill = diss_time), 
                           position = position_dodge(width = 0.8)) + ylim(0,ylim) + 

                geom_text(data = subset(pvals_4N, diss_time == time), aes(x = title, y = ylim, label = significance, group = diss_time),
                        position = position_dodge(width = 0.8), vjust = 0.5, size = 6) +
        
                scale_color_manual(values = c("E7.5" = stats$facet_col, "E8.5" = stats$facet_col2)) + 
                scale_fill_manual(values = c("E7.5" = stats$facet_col, "E8.5" = stats$facet_col2)) + 
                xlab("") + ylab("Lineage Fraction") +
                theme_minimal(base_size = 18) + theme(plot.title = element_text(size = 25, hjust = 0.5, face = "bold"),
                                                      axis.text.x = element_text(size = 12),
                                                      legend.position = "None")
        show_fig(sprintf("./figs/Fig1/violins/%s_%s_4N.pdf", lin, gsub("//.", "_", time)), p, width = 4, height = 3)
    }
}

produce_2N_violins <- function(lin, embryo_orders, ct_md = cell_type_md){
    
    stats = organize_lineage_stats(lin)
    
    plot_df = stats$plot_df
    error_bars = stats$error_bars
    
    ylim = max(stats$plot_df$total_f) + 0.05
    
    plot_df_2N = plot_df %>% filter(title %in% c("2N TKO", "2N DKO", "TKO host", "DKO host"))
    error_bars_2N = error_bars %>% filter(title %in% c("2N TKO", "2N DKO", "TKO host", "DKO host"))
    
    pvals_2N <- plot_df_2N %>%
        group_by(title, diss_time) %>%
            mutate(
                p_value = case_when(grepl("TKO", title) ~ wilcox.test(total_f, plot_df$total_f[plot_df$title == "TKO host" & plot_df$diss_time == diss_time])$p.value,
                                    grepl("DKO", title) ~ wilcox.test(total_f, plot_df$total_f[plot_df$title == "DKO host" & plot_df$diss_time == diss_time])$p.value
                                    )
            )

        # Step 2: Add significance levels based on p-values
        pvals_2N <- pvals_2N %>%
            mutate(significance = case_when(
                p_value < 0.001 ~ "***",
                p_value < 0.01 ~ "**",
                p_value < 0.05 ~ "*",
                TRUE ~ "ns"
                )) %>%
            filter(!grepl("host", title)) %>% select(title, diss_time, p_value, significance) %>% distinct()
    
        for (time in c("E7.5", "E8.5")){
            p <- ggplot() +
                    geom_errorbar(data = subset(error_bars_2N, diss_time == time), 
                                  aes(x = title, ymin=sd_neg, ymax=sd_pos, color = diss_time), 
                                  width = 0.3, position = position_dodge(width = 0.8)) + 
                    geom_errorbar(data = subset(error_bars_2N, diss_time == time), 
                                  aes(x = title, ymin = median, ymax = median, color = diss_time), 
                                  size = 1.5, width = 0.5, position = position_dodge(width = 0.8)) + 
                    geom_violin(data = subset(plot_df_2N, diss_time == time), 
                                  aes(title, total_f, fill = diss_time), 
                                  position = position_dodge(width = 0.8)) + ylim(0,ylim) + 

                    geom_text(data = subset(pvals_2N, diss_time == time), aes(x = title, y = ylim, label = significance, group = diss_time),
                            position = position_dodge(width = 0.8), vjust = 0.5, size = 6) +

                    scale_fill_manual(values = c("E7.5" = stats$facet_col, "E8.5" = stats$facet_col2)) + 
                    scale_color_manual(values = c("E7.5" = stats$stats$facet_col, "E8.5" = stats$facet_col2)) + 
                    xlab("") + ylab("Lineage Fraction") +
                    theme_minimal(base_size = 18) + theme(plot.title = element_text(size = 25, hjust = 0.5, face = "bold"),
                                                          axis.text.x = element_text(size = 12), legend.position = "None")
            show_fig(sprintf("./figs/Fig1/violins/%s_%s_2N.pdf", lin, gsub("//.", "_", time)), p, width = 4, height = 3)

    }
}

# Loop through the lineages: Ectoderm, Endoderm, Embryonic mesoderm, ExE mesoderm - and produce violin plots for 4N and 2N embryos.
plot_all_violins <- function(lineages = c("Ectoderm", "Endoderm", "Embryonic mesoderm", "ExE mesoderm"), embryo_orders = embryo_orders){
    for (lin in lineages){
        produce_4N_violins(lin, embryo_orders = embryo_orders, ct_md = cell_type_md)
        produce_2N_violins(lin, embryo_orders = embryo_orders, ct_md = cell_type_md)
    }
}

produce_supp_violins <- function(lin, embryo_orders = embryo_orders, ct_md = cell_type_md){
    
    stats = organize_lineage_stats(lin)
    
    plot_df = stats$plot_df
    error_bars = stats$error_bars
    
    ylim = max(stats$plot_df$total_f) + 0.05

    plot_df_2N = plot_df %>% filter(title %in% c("DNMT1 host", "DNMT1 KO", "UHRF1 host", "UHRF1 KO", "DNMT3A host", "DNMT3A KO", "DNMT3B host", "DNMT3B KO")) %>%
        mutate(title = factor(title, c("DNMT1 host", "DNMT1 KO", "UHRF1 host", "UHRF1 KO", "DNMT3A host", "DNMT3A KO", "DNMT3B host", "DNMT3B KO")))
    error_bars_2N = error_bars %>% filter(title %in% c("DNMT1 host", "DNMT1 KO", "UHRF1 host", "UHRF1 KO", "DNMT3A host", "DNMT3A KO", "DNMT3B host", "DNMT3B KO"))
    
    pvals_2N <- plot_df_2N %>%
        group_by(title) %>%
            mutate(
                p_value = case_when(grepl("DNMT1", title) ~ wilcox.test(total_f, plot_df$total_f[plot_df$title == "DNMT1 host"])$p.value,
                                    grepl("UHRF1", title) ~ wilcox.test(total_f, plot_df$total_f[plot_df$title == "UHRF1 host"])$p.value,
                                    grepl("DNMT3A", title) ~ wilcox.test(total_f, plot_df$total_f[plot_df$title == "DNMT3A host"])$p.value,
                                    grepl("DNMT3B", title) ~ wilcox.test(total_f, plot_df$total_f[plot_df$title == "DNMT3B host"])$p.value
                                    )
            )

    # Add significance levels based on p-values
    pvals_2N <- pvals_2N %>%
        mutate(significance = case_when(
            p_value < 0.001 ~ "***",
            p_value < 0.01 ~ "**",
            p_value < 0.05 ~ "*",
            TRUE ~ "ns"
            )) %>%
        filter(!grepl("host", title)) %>% select(title, p_value, significance) %>% distinct()

        p <- ggplot() +
                geom_errorbar(data = subset(error_bars_2N), 
                              aes(x = title, ymin=sd_neg, ymax=sd_pos), color = stats$facet_col, 
                              width = 0.3, position = position_dodge(width = 0.8)) + 
                geom_errorbar(data = subset(error_bars_2N), 
                              aes(x = title, ymin = median, ymax = median), color = stats$facet_col, 
                              size = 1.5, width = 0.5, position = position_dodge(width = 0.8)) + 
                geom_violin(data = subset(plot_df_2N), 
                              aes(title, total_f), fill = stats$facet_col, 
                              position = position_dodge(width = 0.8)) + ylim(0,ylim) + 

                geom_text(data = subset(pvals_2N), aes(x = title, y = ylim, label = significance, group = "E7.5"),
                        position = position_dodge(width = 0.8), vjust = 0.5, size = 6) +

                xlab("") + ylab("Lineage Fraction") +
                theme_minimal(base_size = 18) + theme(plot.title = element_text(size = 25, hjust = 0.5, face = "bold"),
                                                      axis.text.x = element_text(size = 12), legend.position = "None")
        # The original run took the file name's stage from a leftover global `time` ("E8.5"); the plot itself pools all stages.
        show_fig(sprintf("./figs/Supp/violins/%s_E8.5_singles.pdf", lin), p, width = 10, height = 3)

        return(p)
    
}

# Loop through the lineages: Ectoderm, Endoderm, Embryonic mesoderm, ExE mesoderm - and produce Supp violin plots for single mutant embryos.
plot_supp_violins <- function(lineages = c("Ectoderm", "Endoderm", "Embryonic mesoderm", "ExE mesoderm"), embryo_orders = embryo_orders){
    for (lin in lineages){
        produce_supp_violins(lin, embryo_orders = embryo_orders, ct_md = cell_type_md)
    }
}
