# KO vs host pseudo-bulk expression per genotype x lineage of the chimeras.
# Shared by Fig. 1i (Fig1.ipynb), Fig. 2e (Fig2.ipynb), Extended Data Fig. 5a (ExtData_RNA.ipynb)
# and Extended Data Fig. 8c (Fig5.ipynb). The chi-square screen (Supplementary Table 1) is run in
# prep/ko_vs_host.ipynb.
#
# Usage (from the repo root, after source("paths.R")):
#   source(file.path(repo_dir, "scripts/single_cell_func.r"))
#   sc <- prep_single_cell()

library(dplyr)
library(tidyr)
library(tibble)
library(tgstat)
library(tgutil)
library(tglkmeans)  # downsample_matrix()

# Embryo-matched KO and host pseudo-bulk UMIs per genotype x lineage (columns "<genotype>.<lineage>").
# For each genotype and lineage, every embryo with more than 2 KO and more than 2 host cells in that
# lineage contributes its summed KO and host UMIs, downsampled to the smaller of the two totals
# (seed 19); the embryos are then summed into one KO column (mko) and one host column (mhost).
# Cached in data/output/mko.RDS and mhost.RDS; recomputed from the single cells only when missing.
load_ko_host_pseudobulk <- function(recompute = FALSE) {
    if (!recompute && file.exists("./data/output/mko.RDS") && file.exists("./data/output/mhost.RDS")) {
        return(list(mko = readRDS(in_path("./data/output/mko.RDS")), mhost = readRDS(in_path("./data/output/mhost.RDS"))))
    }

    library(anndata)
    message("will read cells")
    cells = read_h5ad("./data/MC_models/RNA_all.h5ad")
    message("done read cells")
    umis = t(as.matrix(cells$X))
    message("done cells 2 mat")

    cells$obs$cell_name = rownames(cells$obs)
    colnames(umis) = cells$obs$cell_name

    cell_host = cells$obs$KO_host
    names(cell_host) = cells$obs$cell_name
    cell_geno = as.character(cells$obs$genotype)
    names(cell_geno) = cells$obs$cell_name
    cell_embko = paste(cell_host, cells$obs$embryo, sep="_")
    names(cell_embko) = cells$obs$cell_name
    cell_emb = as.character(cells$obs$embryo)
    names(cell_emb) = cells$obs$cell_name

    # Each cell gets the cell type of its metacell, and each cell type a broad lineage ("layer")
    cell_annots = read.csv("./data/ref/cell_annotation.csv", stringsAsFactors=F)
    type_layer = read.table("./data/ref/annot_key.txt", stringsAsFactors=F,h=T, sep="\t")

    cell_mc_type = cell_annots$cell_type
    names(cell_mc_type) = cell_annots$cell

    type_to_layer = type_layer$layer
    names(type_to_layer) = type_layer$type

    cell_layer = type_to_layer[cell_mc_type]
    names(cell_layer) = cell_annots$cell

    layers = unique(cell_layer)
    layers = setdiff(layers, "NA")
    layers = layers[!is.na(layers)]
    genotypes = unique(cell_geno[names(cell_layer)])
    genotypes = genotypes[!is.na(genotypes)]

    margs_ko = list()
    margs_host = list()
    for(geno in genotypes) {
        marg_ko = NULL
        marg_host = NULL
        for(layer in layers) {
            fnms = intersect(colnames(umis), names(which(cell_layer==layer)))
            fnms = intersect(fnms, names(which(cell_geno==geno)))
            embs = unique(cell_emb[fnms])

            ds_umis = NULL
            for(emb in embs) {
                emb_ko = sprintf("KO_%s",emb)
                emb_host = sprintf("host_%s",emb)
                nms_ko = names(which(cell_embko[fnms] == emb_ko))
                nms_host = names(which(cell_embko[fnms] == emb_host))
                if(length(nms_ko) > 2 & length(nms_host) > 2) {
                    umi_ko = rowSums(umis[,nms_ko])
                    umi_host = rowSums(umis[,nms_host])
                    a = downsample_matrix(cbind(umi_ko, umi_host), min(sum(umi_ko), sum(umi_host)), seed=19)
                    rownames(a) = names(umi_ko)
                    colnames(a) = c(emb_ko, emb_host)
                    if(is.null(a)) {
                        ds_umis = a
                    } else {
                        ds_umis = cbind(ds_umis, a)
                    }
                }
            }

            message("geno ", geno, " layer ", layer, " total ", sum(ds_umis)/2)
            if(!is.null(ds_umis)){
                uko = rowSums(ds_umis[,seq(1,ncol(ds_umis),2)])
                uhost = rowSums(ds_umis[,seq(2,ncol(ds_umis),2)])
                if(is.null(marg_ko)) {
                    marg_ko = uko
                    marg_host = uhost
                    act_layers = layer
                } else {
                    marg_ko = cbind(marg_ko, uko)
                    marg_host = cbind(marg_host, uhost)
                    act_layers = c(act_layers, layer)
                }
            }
        }
        colnames(marg_ko) = paste(rep(geno,t=length(act_layers)), act_layers, sep=".")
        colnames(marg_host) = colnames(marg_ko)
        margs_ko[[geno]] = marg_ko
        margs_host[[geno]] = marg_host
    }
    mko = do.call('cbind', margs_ko)
    mhost = do.call('cbind', margs_host)

    saveRDS(mko, out_path("./data/output/mko.RDS"))
    saveRDS(mhost, out_path("./data/output/mhost.RDS"))
    list(mko = mko, mhost = mhost)
}

# Strip "_" and "2N." from genotype x lineage names for the Fig. 1i heatmap
fix_names = function(names){
    names = gsub("_", " ", names)
    names = gsub("2N.", "", names)
    return(names)
}

# Everything the single-cell level panels need, as one list:
#   mko, mhost            embryo-matched KO / host pseudo-bulk UMIs (genotype x lineage columns)
#   mko_lf, mhost_lf      Fig. 1i: each column downsampled to 250,000 UMIs, log2(2e-5 + frequency), only
#                         columns with >= 250,000 UMIs on both sides, minus each gene's median across columns
#   mko_fe, mhost_fe      Fig. 2e, ED 5a: frequencies of mko / mhost (not downsampled)
#   mdlt_full             log2(1e-5 + KO freq) - log2(1e-5 + host freq), all genes (Fig. 2e, ED 8c)
#   f_dlt, plotted_genes  genes whose change varies by more than 2 across all columns (plotted in Fig. 2e, ED 5a)
#   weber_germ            germline (upregulated_germline) genes of Dahlet et al., 2020
#   weber_highlight       weber_germ genes present in mdlt_full (blue in Fig. 2e, ED 5a)
#   chisq                 chi-square KO vs host screen (data/output/chisq_screen.RDS)
#   sig_consistent        gene x lineage pairs with q < 0.01 in both Dnmt-TKO and Dnmt-DKO (black in ED 5a)
#   cell_type_ref         lineage colours, with the lineage names used here ("layer")
#   colors                heatmap palette of Fig. 1i
prep_single_cell <- function(recompute_pseudobulk = FALSE) {
    sc = load_ko_host_pseudobulk(recompute = recompute_pseudobulk)
    mko = sc$mko
    mhost = sc$mhost

    # Fig. 1i profiles. Seeds are the ones downsample_matrix() drew in the original session
    # (embdnmt SC_level_analysis.ipynb printed "No seed provided. Using 8056" / "4958");
    # with them Fig. 1i matches the published heatmap cell for cell.
    mko_ds = downsample_matrix(mko, 250000, seed = 8056)
    rownames(mko_ds) = rownames(mko)
    colnames(mko_ds) = colnames(mko)
    mhost_ds = downsample_matrix(mhost, 250000, seed = 4958)
    rownames(mhost_ds) = rownames(mhost)
    colnames(mhost_ds) = colnames(mhost)

    mko_e = t(t(mko_ds)/colSums(mko_ds))
    mhost_e = t(t(mhost_ds)/colSums(mhost_ds))

    mko_l = log2(2e-5+mko_e)
    mhost_l = log2(2e-5+mhost_e)

    genos_to_keep = names(which(colSums(mko) >= 2.5e+5 & colSums(mhost) >= 2.5e+5))
    mko_l = mko_l[,colnames(mko_l) %in% genos_to_keep]
    mhost_l = mhost_l[,colnames(mhost_l) %in% genos_to_keep]

    sc$mko_lf = mko_l - apply(mko_l,1,median)
    sc$mhost_lf = mhost_l - apply(mhost_l,1,median)

    # KO - host changes (Fig. 2e, ED 5a, ED 8c)
    sc$mko_fe = t(t(mko)/colSums(mko))
    sc$mhost_fe = t(t(mhost)/colSums(mhost))
    sc$mdlt_full = log2(1e-5+sc$mko_fe)-log2(1e-5+sc$mhost_fe)
    sc$f_dlt = apply(sc$mdlt_full,1,max)-apply(sc$mdlt_full,1,min)>2
    sc$plotted_genes = names(which(sc$f_dlt))

    # Germline genes of Dahlet et al., 2020 (the rows of mko are the genes of the single-cell matrix)
    tab_weber = read.csv("./data/download/weber_genes.csv",stringsAsFactors=F)
    sc$weber_germ = intersect(rownames(mko), tab_weber$gene[tab_weber$class=="upregulated_germline"])
    sc$weber_highlight = intersect(rownames(sc$mdlt_full), sc$weber_germ)

    # Chi-square screen (made in prep/ko_vs_host.ipynb)
    sc$chisq = readRDS(in_path("./data/output/chisq_screen.RDS"))
    sc$sig_consistent = sc$chisq %>% group_by(gene, layer) %>% mutate(sum_sig = sum(qval < 0.01)) %>% filter(sum_sig == 2) %>%
        select(gene, layer) %>% distinct()

    # Lineage colours
    cell_type_ref = read.csv("./data/ref/cell_type_ref.csv")
    cell_type_ref = cell_type_ref %>% select(lineage, lineage_color) %>% distinct() %>% filter(
        lineage %in% c("Endoderm", "Endothelial cells", "Blood", "ExE mesoderm", "Embryonic mesoderm", "Epiblast", "Ectoderm"))
    cell_type_ref$layer = c("Endoderm", "Blood", "Endothel", "ExMeso", "EmbMeso", "Epi", "Ecto")
    sc$cell_type_ref = cell_type_ref

    colorscale = c("black", "grey30", "grey60", "grey90", "white", "deepskyblue1", "deepskyblue2", "deepskyblue3", "deepskyblue4")
    sc$colors = colorRampPalette(colorscale)(100)

    sc
}

# Chi-square KO vs host screen (Supplementary Table 1): for every gene, genotype and lineage, a 2 x 2 test of
# the gene's UMIs and all other UMIs in KO vs host cells; genes with < UMI_thresh UMIs (KO + host) skipped,
# FDR per lineage.
calc_chisq_screen <- function(mko, mhost, genos = c("DNMT_TKO", "DNMT_DKO", "TET_TKO"), UMI_thresh = 20, num_proc = 40) {
    ko_for_chi = mko %>% as.data.frame() %>% rownames_to_column("gene") %>% pivot_longer(-gene, names_to = "geno_layer", values_to = "ko") %>%
        group_by(geno_layer) %>% mutate(tot_min_ko = sum(ko)-ko) %>% ungroup()

    host_for_chi = mhost %>% as.data.frame() %>% rownames_to_column("gene") %>% pivot_longer(-gene, names_to = "geno_layer", values_to = "host") %>%
        group_by(geno_layer) %>% mutate(tot_min_host = sum(host)-host) %>% ungroup()

    chi_df = full_join(ko_for_chi, host_for_chi) %>% mutate(geno = sapply(strsplit(geno_layer, "\\."), `[`, 1),
                                                            layer = sapply(strsplit(geno_layer, "\\."), `[`, 2)) %>%
        select(-geno_layer) %>% filter(geno %in% genos)

    doMC::registerDoMC(num_proc)
    res <- chi_df %>%
        filter(!(ko + host < UMI_thresh),
               rowSums(. < 0) == 0) %>%
        plyr::adply(1, function(x) {
            chisq <- chisq.test(
                matrix(c(x$ko, x$tot_min_ko, x$host, x$tot_min_host), byrow = TRUE, nrow = 2)
            )
            tibble(pval = chisq$p.value, statistic = chisq$statistic)
        }, .parallel = TRUE) %>%
        group_by(layer) %>%
        mutate(qval = p.adjust(pval, method = "fdr")) %>%
        ungroup()

    res %>% mutate(dir = case_when(ko > host ~ "up",
                                   ko < host ~ "down",
                                   ko == host ~ "same"),
                   sig = qval < 0.01) %>%
        select(gene, geno, layer, statistic, qval, dir, sig)
}
