library(scales)
library(colorspace)
library(misha.ext)
library(colorspace)
library(tidyr)
library(tglkmeans)
library(pheatmap)
library(data.table)
library(mcATAC)
library(zoo)
library(dplyr)
library(ggplot2)
library(tibble)

init_e95 <- function(mod_prev, force_cmat_redo = FALSE){

    mod_e95 <- list()
    base_dir = paste0(embdnmt_dir, "/revision_e95")
    setwd(base_dir)

    if (file.exists(sprintf("%s/data/output/peaks_on_sum.RDS", base_dir))){
        peaks_on_sum <- readRDS(in_path(sprintf("%s/data/output/peaks_on_sum.RDS", base_dir)))
        mod_e95$npeaks <- peaks_on_sum
    } else {
        peaks_on_sum = comp_peaks_over_punc_region("embdnmt.brain_e95_dko.marginal_plus_gastr",
                                                   "embdnmt.brain_e95_dko.marginal_plus_gastr_regcap",
                                                   peak_ext = 140, eps = 3, T_q = 0.97)
        rownames(peaks_on_sum) = peak_names(peaks_on_sum)
        peaks_on_sum = peaks_on_sum %>% filter(chrom != "chr14")
        saveRDS(peaks_on_sum, out_path(sprintf("%s/data/output/peaks_on_sum.RDS", base_dir)))
        mod_e95$npeaks <- peaks_on_sum
    }

    # 
    peaks = mod_e95$npeaks
    peaks$peak_name = rownames(peaks)
    tss_peaks = gintervals.neighbors(peaks, mod_prev$tss, mindist = -5e+3, maxdist = 5e+3) %>% select(chrom, start, end, peak_name)
    rownames(tss_peaks) = tss_peaks$peak_name
    cre_peaks = peaks[setdiff(rownames(peaks), rownames(tss_peaks)),]

    # Normalize for energy calculations, make sure no overlapping regions
    cre_peaks_for_energ = gintervals.normalize(cre_peaks, 500)
    union_on_overlapping = gintervals.union(cre_peaks_for_energ, cre_peaks_for_energ)
    cre_peaks_for_energ = gintervals.normalize(union_on_overlapping, 500)
    cre_peaks_for_energ = gintervals.neighbors(cre_peaks_for_energ, cre_peaks, mindist = 0, maxdist = 0, maxneighbors = 1) %>% select(chrom, start, end, peak_name)
    mod_e95$cre_peaks_for_energ = cre_peaks_for_energ

    tss_peaks_for_energ = gintervals.normalize(tss_peaks, 500)
    union_on_overlapping = gintervals.union(tss_peaks_for_energ, tss_peaks_for_energ)
    tss_peaks_for_energ = gintervals.normalize(union_on_overlapping, 500)
    tss_peaks_for_energ = gintervals.neighbors(tss_peaks_for_energ, tss_peaks, mindist = 0, maxdist = 0, maxneighbors = 1) %>% select(chrom, start, end, peak_name)

    cre_peaks = cre_peaks[cre_peaks_for_energ$peak_name,] %>% distinct() 
    tss_peaks = tss_peaks[tss_peaks_for_energ$peak_name,] %>% distinct()
    
    mod_e95$tss_peaks = tss_peaks
    mod_e95$cre_peaks = cre_peaks

    # Load all RDSs, if they exist, otherwise create them
    if (file.exists(sprintf("%s/data/output/cov_chr.RDS", base_dir))){
        cov_chr <- readRDS(in_path(sprintf("%s/data/output/cov_chr.RDS", base_dir)))
        mod_e95$cov_chr <- cov_chr
    } else {
        collect_cov_r.e95(mod_e95$npeaks)
        cov_chr <- readRDS(in_path(sprintf("%s/data/output/cov_chr.RDS", base_dir)))
        mod_e95$cov_chr <- cov_chr
    }

    states_host = cov_chr[[1]][,1]
    states_host = states_host[!is.na(states_host)]
    states_dko = cov_chr[[2]][,1]
    states_dko = states_dko[!is.na(states_dko)]

    states_both = intersect(names(states_host), names(states_dko))
    states_both = states_both[states_both != "marginal"]
    mod_e95$states_both <- states_both

    mod_e95$ecto_states = c("Hindbrain", "Optic_vesicle_Telencephalon", "Midbrain_hindbrain_boundary", 
                            "Dmbx1_Diencephalon_Mesencephalon", "Hoxd3_spinal_cord_progenitors", "Neural_crest_1")
    mod_e95$meso_states = c("Somitic", "Sclerotome", "Bmp4_Igf2_Cranial_muscles")

    if (file.exists(sprintf("%s/data/output/multi_type_cov.RDS", base_dir))){
        cov_stat <- readRDS(in_path(sprintf("%s/data/output/multi_type_cov.RDS", base_dir)))
        mod_e95$cov_stat <- cov_stat
    }

    mod_e95$cov_stat_kary_norm <- generate_stat_kary_norm.e95(mod_e95)
    mod_e95$gastr_states = setdiff(mod_prev$cov_states, "marginal")
    mod_e95$brain_ct_order = c('Hindbrain.host', 'Hindbrain.dko', 
                'Hoxd3_spinal_cord_progenitors.host', 'Hoxd3_spinal_cord_progenitors.dko', 
                'Midbrain_hindbrain_boundary.host', 'Midbrain_hindbrain_boundary.dko', 
                'Optic_vesicle_Telencephalon.host', 'Optic_vesicle_Telencephalon.dko', 
                'Dmbx1_Diencephalon_Mesencephalon.host', 'Dmbx1_Diencephalon_Mesencephalon.dko')
    
    mod_e95$cmat_nc <- generate_cmat(mod_e95, mod_prev, force_redo = force_cmat_redo)
    mod_e95$cmat <- generate_cmat(mod_e95, mod_prev, 
                                  states_to_exclude = c("Neural_crest_1"), 
                                  force_redo = force_cmat_redo)

    gastr_state_key = read.csv(paste0(embdnmt_dir, "/data/ref/cell_type_key.csv"))
    mod_e95$gastr_state_colors = gastr_state_key %>% select(shorthand, color) %>% deframe()

    mod_e95$energ_cre_fulldb = readRDS(in_path(sprintf("%s/data/output/energ_cre_fulldb.RDS", base_dir)))
    mod_e95$energ_cre_blm = readRDS(in_path(sprintf("%s/data/output/energ_cre_blm.RDS", base_dir)))
    
    colors_e95 = c("Hindbrain" = "#4a7c59", 
               "Hoxd3+ (spinal cord progenitors)" = "#68b0ab", 
               "Dorsal midbrain" = "#8fc0a9", 
               "Midbrain–hindbrain boundary" = "#216869", 
               "Optic vesicle / Telencephalon" = "#faf3dd", 
               "Dmbx1+ (Diencephalon/Mesencephalon)" = "#c8d5b9",
              
               "Neural crest 1" = "#93da49",
               "Neural crest 2" = "#c2e812",
               
               "Epidermis" = "#ffd6ba",
               "Endothelium" = "#ff891c",
               "Somitic" = "#62b6cb",
               "Sclerotome" = "#90d0da",
               "Bmp4+ Igf2+ / Cranial muscles" = "#c6c4f0",
               "Cardiomyocytes" = "#b51d8d",
               "Erythroid" = "#c72228",
               "Megakaryocytes" = "#FF00FF")
    
    colors_e95_flat_names = setNames(unname(colors_e95), convert_names(names(colors_e95)))
    
    mod_e95$colors_e95 <- colors_e95
    mod_e95$colors_e95_flat_names <- colors_e95_flat_names
    
    return(mod_e95)
}

convert_names <- function(vec){
        vec = gsub(" ", "_", vec)
        vec = gsub("\\/", "_", vec)
        vec = gsub("\\+", "_", vec)
        vec = gsub("[()]", "", vec)
        vec = gsub("\\p{Pd}+", "_", vec, perl = TRUE)
        vec = gsub("_+", "_", vec)

        return(vec)
}

collect_cov_r.e95 <- function(npeaks) 
{
    tss = gintervals.load("intervs.global.tss")
    tss = tss[tss$chrom != "chrM" & tss$chrom != "chrY", ]
    all_nms = gtrack.ls("embdnmt.brain_e95_host.reannot.")
    all_nms = union(all_nms, gtrack.ls("embdnmt.brain_e95_dko.reannot."))
    all_nms = sub("embdnmt.brain_e95_host.reannot.", "", all_nms)
    all_nms = sub("embdnmt.brain_e95_dko.reannot.", "", all_nms)
    all_nms = all_nms[!grepl("_normed", all_nms) & !grepl("_smoothed", all_nms)]
    all_nms = unique(all_nms)
    
    f_rp = grepl("^Rp[l|s]", tss$geneSymbol)
    rp_tss = tss[f_rp, c("chrom", "start", "end", "strand")]
    rp_tss$start = rp_tss$start + ifelse(rp_tss$strand == 1, 
        -800, -200)
    rp_tss$end = rp_tss$end + ifelse(rp_tss$strand == 1, 200, 
        800)
    cov_stats = list()
    cov_peak_chr = list()
    cov_chr = list()
    for (ds in c("brain_e95_host.reannot", "brain_e95_dko.reannot")) {
        cov_stat = NULL
        for (nm in all_nms) {
            message("doing track ", ds, " ", nm)
            track_nm = sprintf("embdnmt.%s.%s", ds, nm)
            if (gtrack.exists(track_nm)) {
         
                tot = as.numeric(gsummary(track_nm)[["Sum"]])
                a = gextract(track_nm, intervals = npeaks, iterator = 20, 
                  colnames = c("cov"))
                a_rp = gextract(track_nm, intervals = rp_tss, 
                  iterator = 20, colnames = c("cov"))
                gvtrack.create("tcov", track_nm, "sum")
                a_chrom = gextract(c("tcov"), intervals = gintervals.all(), 
                  iterator = gintervals.all(), colnames = c("cov"))
                tot_peak = sum(a$cov, na.rm = T)
                tot_peak_chr = tapply(a$cov, a$chrom, sum, na.rm = T)
                tot_rp = sum(a_rp$cov, na.rm = T)
                
                if (is.null(cov_stat)) {
                  cov_stat = data.frame(nm = nm, tot = tot, tot_peak = tot_peak, 
                    tot_rp = tot_rp)
                  cov_peak_chr_st = c(nm, tot_peak_chr)
                  cov_chr_st = c(nm, a_chrom$cov)
                }
                else {
                  cov_stat = rbind(cov_stat, c(nm, tot, tot_peak, 
                    tot_rp))
                  cov_peak_chr_st = rbind(cov_peak_chr_st, c(nm, 
                    tot_peak_chr))
                  cov_chr_st = rbind(cov_chr_st, c(nm, a_chrom$cov))
                }
            }
            else {
                if (is.null(cov_stat)) {
                  cov_stat = data.frame(nm = nm, tot = NA, tot_peak = NA, 
                    tot_rp = NA)
                  cov_peak_chr_st = c(nm, rep(NA, 22))
                  cov_chr_st = c(nm, rep(NA, 22))
                }
                else {
                  cov_stat = rbind(cov_stat, c(nm, NA, NA, NA))
                  cov_peak_chr_st = rbind(cov_peak_chr_st, c(nm, 
                    rep(NA, 22)))
                  cov_chr_st = rbind(cov_chr_st, c(nm, rep(NA, 
                    22)))
                }
            }
        }
        cov_stat$tot = as.numeric(cov_stat$tot)
        cov_stat$tot_peak = as.numeric(cov_stat$tot_peak)
        colnames(cov_stat)[2] = sprintf("%s.tot", ds)
        colnames(cov_stat)[3] = sprintf("%s.tot_peak", ds)
        colnames(cov_stat)[4] = sprintf("%s.tot_rp", ds)
        cov_stats[[ds]] = cov_stat
        rownames(cov_peak_chr_st) = cov_peak_chr_st[, 1]
        cov_peak_chr_st = cov_peak_chr_st[, -1]
        class(cov_peak_chr_st) = "numeric"
        cov_peak_chr[[ds]] = cov_peak_chr_st
        rownames(cov_chr_st) = cov_chr_st[, 1]
        cov_chr_st = cov_chr_st[, -1]
        colnames(cov_chr_st) = gintervals.all()$chrom
        class(cov_chr_st) = "numeric"
        cov_chr[[ds]] = cov_chr_st
    }
    cov_stat = cov_stats[["brain_e95_host.reannot"]] %>% left_join(cov_stats[["brain_e95_dko.reannot"]])

    base_dir = paste0(embdnmt_dir, "/revision_e95")
    saveRDS(cov_stat, file = out_path(sprintf("%s/data/output/multi_type_cov.RDS", base_dir)))
    saveRDS(cov_peak_chr, file = out_path(sprintf("%s/data/output/cov_peak_chr.RDS", base_dir)))
    saveRDS(cov_chr, file = out_path(sprintf("%s/data/output/cov_chr.RDS", base_dir)))
    rownames(cov_stat) = cov_stat$nm
}

gen_track_vext.e95 <- function(tname, w_ext = 140){
    message("set vt for ", tname)
    nm_host = sprintf("embdnmt.brain_e95_host.reannot.%s", tname)
    nm_dko = sprintf("embdnmt.brain_e95_dko.reannot.%s", tname)
    gvtrack.create("host_ext", nm_host, "sum")
    gvtrack.create("dko_ext", nm_dko, "sum")
    gvtrack.iterator("host_ext", sshift = -w_ext, eshift = w_ext)
    gvtrack.iterator("dko_ext", sshift = -w_ext, eshift = w_ext)
}

build_peak_mat.e95 = function(npeaks, state, cov_stat_kary_norm)
{
    gen_track_vext.e95(state, w_ext=140)
    
	bg_host = "embdnmt.brain_e95_host.marginal_regcap"
	bg_dko = "embdnmt.brain_e95_dko.marginal_regcap"
	mat = gextract(c("host_ext", "dko_ext", bg_host, bg_dko),
				   intervals = npeaks, iterator = npeaks, 
				   colnames = c("host_n", "dko_n", "bg_host", "bg_dko"))

	mat[is.na(mat)] = 0

    
    
    # decide on ds coverage (these numbers already account for karyo correction)
    r_host = as.numeric(cov_stat_kary_norm[state,"brain_e95_host.reannot.tot_rp"])
	r_dko = as.numeric(cov_stat_kary_norm[state,"brain_e95_dko.reannot.tot_rp"])
    r_min = min(r_host, r_dko)
    
	r_host = r_min/r_host
	r_dko = r_min/r_dko
    
    # correct standartized coverage
	if(r_host < 1) {
		to_ds = matrix(mat$host_n,ncol=1)
		post_ds = downsample_matrix(to_ds, target_n = floor(sum(to_ds)*r_host), seed=19)
		mat$host_ds = post_ds
	} else {
		mat$host_ds = mat$host_n
	}
	if(r_dko < 1) {
		to_ds = matrix(mat$dko_n,ncol=1)
		post_ds = downsample_matrix(to_ds, target_n = floor(sum(to_ds)*r_dko), seed=19)
		mat$dko_ds = post_ds
	} else {
		mat$dko_ds = mat$dko_n
	}

	mat$host = mat$host_n * r_host
	mat$dko = mat$dko_n * r_dko

	return(mat)	
}

generate_stat_kary_norm.e95 = function(mod_e95){
    cov_chr = mod_e95$cov_chr
    cov_stat = mod_e95$cov_stat

    r_chr_cov_dko = cov_chr[["brain_e95_dko.reannot"]]/cov_chr[["brain_e95_host.reannot"]]
    r_chr_cov_dko = r_chr_cov_dko/apply(r_chr_cov_dko, 1, median, na.rm=T)

    dko_tot_karyo = rowSums(cov_chr[["brain_e95_dko.reannot"]],na.rm=T)/rowSums(cov_chr[["brain_e95_dko.reannot"]]/r_chr_cov_dko, na.rm=T)
        
    cov_norm = cov_stat
    cov_norm$dko.tot = floor(cov_stat$brain_e95_dko.reannot.tot / dko_tot_karyo)
    cov_norm$dko.tot_peak = floor(cov_stat$brain_e95_dko.reannot.tot_peak / dko_tot_karyo)
    cov_norm$dko.tot_rp = floor(as.numeric(cov_stat$brain_e95_dko.reannot.tot_rp) / dko_tot_karyo)

    cov_stat_kary_norm = cov_norm
    rownames(cov_stat_kary_norm) = cov_stat_kary_norm$nm

    mod_e95$cov_stat_kary_norm = cov_stat_kary_norm
    return(mod_e95)
}

generate_cmat <- function(mod_e95, mod_prev, states_to_exclude = NULL, force_redo = F, reg = 1e-6){

    cov_chr <- mod_e95$cov_chr
    states_to_crecl = c(mod_e95$ecto_states, mod_e95$meso_states)
    if (!is.null(states_to_exclude)){
        states_to_crecl = setdiff(states_to_crecl, states_to_exclude)
    }
    npeaks = mod_e95$npeaks
    cov_stat_kary_norm = mod_e95$cov_stat_kary_norm$cov_stat_kary_norm

    gastr_states = mod_e95$gastr_states
    gastr_state_colors = mod_e95$gastr_state_colors

    base_dir = paste0(embdnmt_dir, "/revision_e95")
    if (file.exists(sprintf("%s/data/output/gastr_crecl_mat.RDS", base_dir))){
        gastr_crecl_mat = readRDS(in_path(sprintf("%s/data/output/gastr_crecl_mat.RDS", base_dir)))
    } else {
        crecl_mat = NULL
        for (st in gastr_states){

            message("doing state ", st)

            pmat = build_peak_mat(mod = mod_prev, peaks = npeaks, state = st)

            pmat$host_p = pmat$wt_n/sum(pmat$wt_n)
            pmat$dko_p = pmat$dnmtd_n/sum(pmat$dnmtd_n)

            pmat$host_p = log2(pmat$host_p+reg)
            pmat$dko_p = log2(pmat$dko_p+reg)

            insert = pmat[,c("host_p", "dko_p")]
            colnames(insert) = c(paste0(st, ".host"), paste0(st, ".dko"))

            if (is.null(crecl_mat)){
                crecl_mat = insert
            } else {
                crecl_mat = cbind(crecl_mat, insert)
            }
        }
        
        saveRDS(crecl_mat, out_path(sprintf("%s/data/output/gastr_crecl_mat.RDS", base_dir)))
    }

    if (file.exists(sprintf("%s/data/output/crecl_mat.RDS", base_dir)) & force_redo == F){
        crecl_mat = readRDS(in_path(sprintf("%s/data/output/crecl_mat.RDS", base_dir)))
    } else {
        crecl_mat = NULL
        for (st in states_to_crecl){

            message("doing state ", st)

            pmat = build_peak_mat.e95(npeaks, st, cov_stat_kary_norm)

            pmat$host_p = pmat$host_n/sum(pmat$host_n)
            pmat$dko_p = pmat$dko_n/sum(pmat$dko_n)

            pmat$host_p = log2(pmat$host_p+reg)
            pmat$dko_p = log2(pmat$dko_p+reg)

            insert = pmat[,c("host_p", "dko_p")]
            colnames(insert) = c(paste0(st, ".host"), paste0(st, ".dko"))

            if (is.null(crecl_mat)){
                crecl_mat = insert
            } else {
                crecl_mat = cbind(crecl_mat, insert)
            }

        }
        
        saveRDS(crecl_mat, out_path(sprintf("%s/data/output/crecl_mat.RDS", base_dir)))
    }

    gastr_crecl_mat = readRDS(in_path(sprintf("%s/data/output/gastr_crecl_mat.RDS", base_dir)))
    crecl_mat = readRDS(in_path(sprintf("%s/data/output/crecl_mat.RDS", base_dir)))

    # Remove states with low coverage
    gastr_crecl_mat = gastr_crecl_mat[,setdiff(colnames(gastr_crecl_mat), c("ac.host", "ac.dko", "pmeso.host", "pmeso.dko"))]
    mod_e95$gastr_crecl_mat = gastr_crecl_mat

    cmat = cbind(crecl_mat, gastr_crecl_mat)
    # Reorder columns
    if (!is.null(states_to_exclude)){
        cell_type_order = c(setdiff(mod_e95$brain_ct_order, states_to_exclude), colnames(crecl_mat)[11:ncol(crecl_mat)], colnames(gastr_crecl_mat))
    } else {
        cell_type_order = c(mod_e95$brain_ct_order, colnames(crecl_mat)[11:ncol(crecl_mat)], colnames(gastr_crecl_mat))
    }
    
    cmat = cmat[,cell_type_order]
    rownames(cmat) = rownames(npeaks)

    return(cmat)
}

filter_peaks_by_diff <- function(mat, diff_thresh, nonactive_thresh = NULL) {
    
    idx <- seq(2, ncol(mat), by = 2)
    diff_mat <- do.call(cbind, lapply(idx, function(i) mat[, i] - mat[, i - 1]))
    rownames(diff_mat) <- rownames(mat)

    keep_diff <- apply(abs(diff_mat) >= diff_thresh, 1, any)

    if (!is.null(nonactive_thresh)) {
    keep_active <- apply(mat, 1, max, na.rm = TRUE) >= nonactive_thresh
    } else {
    keep_active <- rep(TRUE, nrow(mat))
    }

    message("keep_diff: ", sum(keep_diff),
          " | keep_active: ", sum(keep_active),
          " | final: ", sum(keep_diff & keep_active))

    mat_f = mat[keep_diff & keep_active, , drop = FALSE]
    
    return(mat_f)
}

plot_crecl_htmp <- function(mod_e95, cmat_name, peak_subset, n_clust, lower, upper, diff_thresh, passed_thresh = TRUE, nonactive_thresh = NULL, gaps_col = NULL, name){

    cmat = mod_e95[[cmat_name]]
    message("Plotting crecl heatmap for ", name)
    peaks_to_plot = mod_e95[[peak_subset]]
    htmp_mat = cmat[peaks_to_plot$peak_name,]
    message(paste0("n = ", nrow(htmp_mat)))

    htmp_mat = filter_peaks_by_diff(htmp_mat, diff_thresh, nonactive_thresh)
    message(paste0("post threshold: n = ", nrow(htmp_mat)))

    # ANNOTATION and COLORS #
    ann_cols <- list(
    source = c(host = "grey50", dko = "limegreen")
    )
    
    ann_cols$cell_type = mod_e95$colors_e95
    ann_cols$cell_type = c(ann_cols$cell_type, mod_e95$gastr_state_colors)
    names(ann_cols$cell_type) = convert_names(names(ann_cols$cell_type))

    src <- ifelse(grepl("\\.host$", colnames(htmp_mat)), "host",
        ifelse(grepl("\\.dko$",  colnames(htmp_mat)), "dko", NA))

    md = read.csv(in_path("./data/md_metacells.csv"), row.names = 1)
    annotation_dict = deframe(distinct(md[,c("cell_type", "reannot")]))
    
    ann_df <- data.frame(source = src, cell_type = gsub(".host|.dko", "", colnames(cmat)), row.names = colnames(cmat))
    
    cols <- c("white","azure1","azure2","deepskyblue1","deepskyblue3","#0072B2","#08306B","black")
    my_pal <- colorRampPalette(cols)(100)

    # CLUSTERING #
    km = TGL_kmeans_tidy(htmp_mat, k = n_clust, hclust_intra_clusters = TRUE, seed = 30)
    cluster <- km$cluster
    cluster = deframe(cluster)

    centers <- km$centers
    clu_ord <- hclust(dist(centers))$order
    cluster_levels <- as.integer(clu_ord)

    # Row order based on intra cluster ordering
    row_order = c()
    for (cl in cluster_levels){
        cl_peaks_ordered = km$order %>% filter(clust == cl) %>% arrange(intra_clust_order) %>% pull(id)
        row_order = c(row_order, cl_peaks_ordered)
    }

    # PLOTTING #
    htmp_ord <- htmp_mat[row_order, , drop = FALSE]

    ann_rows <- data.frame(
        cluster = factor(cluster[row_order], levels = cluster_levels)
    )

    gaps <- cumsum(table(ann_rows$cluster))

    nbins    <- length(my_pal)
    brks     <- seq(lower, upper, length.out = nbins + 1)
   
    # Cluster run boundaries
    runs <- rle(as.integer(ann_rows$cluster))$lengths
    starts <- c(1, head(cumsum(runs) + 1, -1))
    ends   <- cumsum(runs)
    mids   <- floor((starts + ends) / 2)

    row_lab <- rep("", nrow(htmp_ord))
    clus_ids <- as.integer(levels(ann_rows$cluster))[ann_rows$cluster[mids]]
    row_lab[mids] <- paste0("C", clus_ids)

    htmp_smoothed <- rollmean_in_clusters(mat = as.matrix(htmp_ord), cluster, k = 20)
    show_fig(sprintf("./figs/crecl/%s_crecls.png", name), draw = pheatmap::pheatmap(
        htmp_smoothed,
        cluster_rows = FALSE,
        cluster_cols = FALSE,
        gaps_row = cumsum(runs)[-length(runs)],
        gaps_col = gaps_col,
        color = my_pal, breaks = brks,
        annotation_row = NULL,
        annotation_col = ann_df,
        annotation_colors = ann_cols,
        annotation_legend = FALSE, 
        legend = TRUE,
        na_col = "grey90",
        show_rownames = TRUE,
        show_colnames = FALSE,
        labels_row = row_lab,
        fontsize_row = 15,
        border_color = "black",
        width = 15, height = 15,
        fontsize = 15
        ))

    saveRDS(ann_rows, out_path(sprintf("./data/output/%s_cluster_assignments.RDS", name)))
    saveRDS(ann_cols, out_path(sprintf("./data/output/%s_cluster_colors.RDS", name)))
    saveRDS(htmp_ord, out_path(sprintf("./data/output/%s_cluster_mat.RDS", name)))
    saveRDS(htmp_smoothed, out_path(sprintf("./data/output/%s_cluster_mat_smoothed.RDS", name)))

}

plot_crecl_clusters <- function(
    name,
    clusters_keep = NULL,
    keep_cluster_order = FALSE,      # NEW: keep cluster order exactly like clusters_keep
    out_dir = "./data/output",
    figs_dir = "./figs/crecl"
    ) {

    ann_rows <- readRDS(in_path(file.path(out_dir, sprintf("%s_cluster_assignments.RDS", name))))
    ann_cols <- readRDS(in_path(file.path(out_dir, sprintf("%s_cluster_colors.RDS", name))))
    mat      <- readRDS(in_path(file.path(out_dir, sprintf("%s_cluster_mat_smoothed.RDS", name))))

    cl_full <- suppressWarnings(as.integer(as.character(ann_rows$cluster)))
    if (all(is.na(cl_full))) cl_full <- as.integer(ann_rows$cluster)

    if (!is.null(clusters_keep)) {
      cl <- clusters_keep
      if (is.character(cl)) cl <- as.integer(gsub("^C", "", cl))
      keep <- cl_full %in% as.integer(cl)

      mat <- mat[keep, , drop = FALSE]
      cl_full <- cl_full[keep]

      if (isTRUE(keep_cluster_order)) {
        lev <- as.integer(cl)
        ord <- order(factor(cl_full, levels = lev), seq_along(cl_full))
        mat <- mat[ord, , drop = FALSE]
        cl_full <- cl_full[ord]
      }
    }

    # gaps + labels
    runs <- rle(cl_full)$lengths
    gaps_row <- if (length(runs) > 1) cumsum(runs)[-length(runs)] else NULL

    row_lab <- rep("", nrow(mat))
    r <- rle(cl_full)
    starts <- c(1, head(cumsum(r$lengths) + 1, -1))
    ends   <- cumsum(r$lengths)
    mids   <- floor((starts + ends) / 2)
    row_lab[mids] <- paste0("C", r$values)

    # column annotation (same logic as original)
    src <- ifelse(grepl("\\.host$", colnames(mat)), "host",
           ifelse(grepl("\\.dko$", colnames(mat)), "dko", NA))

    md = read.csv(in_path("./data/md_metacells.csv"), row.names = 1)
    annotation_dict = deframe(distinct(md[,c("cell_type", "reannot")]))

    ann_df <- data.frame(source = src, cell_type = gsub(".host|.dko", "", colnames(mat)), row.names = colnames(mat))

    # palette like original
    cols <- c("white","azure1","azure2","deepskyblue1","deepskyblue3","#0072B2","#08306B","black","#240046","#3c096c")
    my_pal <- colorRampPalette(cols)(100)
    brks <- seq(-20, -15, length.out = length(my_pal) + 1)
    
    suffix <- if (is.null(clusters_keep)) "ALL" else paste0("C", unique(cl_full), collapse = "_")
    out_pdf <- file.path(figs_dir, sprintf("%s_clusters_%s.pdf", name, suffix))
    out_colorscale <- file.path(figs_dir, sprintf("%s_colorscale_new.pdf", name))
    
        # Plot color scale
        z <- matrix(seq_along(my_pal), nrow = length(my_pal), ncol = 1)

        show_fig(out_colorscale, show = FALSE, width = 1.6, height = 6, draw = {
        par(mar = c(0, 6, 0, 1))  # tight, vertical

        plot.new()
        plot.window(xlim = c(0, 1), ylim = c(0, 100))

        rect(
          xleft   = 0,
          ybottom = 0:99,
          xright  = 1,
          ytop    = 1:100,
          col     = my_pal,
          border  = NA
        )

        axis(
          side = 2,
          at = seq(0, 100, length.out = 5),
          labels = round(seq(-20,
                             -14,
                             length.out = 5), 2),
          las = 1,
          cex.axis = 1.2
        )
        })
    
    # Reorder columns
    mat <- mat[,c(17,18,29:30,19:28,31:32,1:16)]

    show_fig(out_pdf, draw = pheatmap::pheatmap(
        mat,
        cluster_rows = FALSE,
        cluster_cols = FALSE,
        gaps_row = gaps_row,
        gaps_col = c(2,4,14,16,26),
        color = my_pal, breaks = brks,
        annotation_col = ann_df,
        annotation_colors = ann_cols,
        annotation_legend = FALSE,
        legend = TRUE,
        na_col = "grey90",
        show_rownames = TRUE,
        show_colnames = FALSE,
        labels_row = row_lab,
        fontsize_row = 15,
        border_color = "black",
        width = 15, height = nrow(mat)/1e+3,
        fontsize = 15
    ))

    invisible(out_pdf)
}

rollmean_in_clusters <- function(mat, cluster, k = 20, align = "center", partial = TRUE) {
    cl <- unname(cluster[rownames(mat)])

    runs <- rle(cl)
    out  <- matrix(NA_real_, nrow(mat), ncol(mat), dimnames = dimnames(mat))

    from <- 1L
    for (i in seq_along(runs$lengths)) {
        to <- from + runs$lengths[i] - 1L
        block <- mat[from:to, , drop = FALSE]

        out[from:to,] <- zoo::rollapply(
            data = block,
            width = k,
            by = 1,
            FUN = function(x) colMeans(x),
            by.column = FALSE,
            align = align,
            partial = partial
        )
        from <- to + 1L
    }
    
    return(out)
}

                                      

# ---- used by E9.5_multiome.Fig6/Fig6.ipynb (Fig. 6e-g) ----

plot_clusters_host_vs_ko_4panels <- function(
    cmat,
    peaks_df,
    name,
    clusters_keep = NULL,
    types_all,
    cluster_assignments_rds = "./data/output/%s_cluster_assignments.RDS",
    point_alpha = 1,
    point_size = 0.1,
    peak_names = NULL,
    use_peak_names_only = TRUE){

    ann_rows <- readRDS(in_path(sprintf(cluster_assignments_rds, name)))
    cl_full <- suppressWarnings(as.integer(as.character(ann_rows$cluster)))
    if (all(is.na(cl_full))) cl_full <- as.integer(ann_rows$cluster)

    normalize_clusters <- function(x) {
        if (is.null(x) || length(x) == 0) return(NULL)
        x <- as.character(x)
        x <- x[x != ""]
        if (!length(x)) return(NULL)
        if (all(grepl("^C\\d+$", x))) x <- gsub("^C", "", x)
        suppressWarnings(as.integer(x))
    }

    cl <- normalize_clusters(clusters_keep)

    # If supplied peak names
    if (!is.null(peak_names) && use_peak_names_only) {
        keep_peaks <- intersect(as.character(peak_names), rownames(cmat))
    } else {
            if (is.null(cl) || all(is.na(cl))) {
                keep_peaks <- intersect(rownames(ann_rows), rownames(cmat))
        } else {
            keep_peaks <- intersect(rownames(ann_rows)[cl_full %in% cl], rownames(cmat))
        }
        if (!is.null(peak_names)) {
            keep_peaks <- intersect(keep_peaks, as.character(peak_names))
        }
    }

    host_cols <- vapply(types_all, `[[`, "", "host_col")
    ko_cols   <- vapply(types_all, `[[`, "", "ko_col")

    get_mean_df <- function(idx, label) {
        tibble::tibble(
          peak_name = keep_peaks,
          group = label,
          host = rowMeans(cmat[keep_peaks, host_cols[idx], drop = FALSE], na.rm = TRUE),
          ko   = rowMeans(cmat[keep_peaks, ko_cols[idx],   drop = FALSE], na.rm = TRUE)
        )
    }
    get_single_df <- function(idx, label) {
        tibble::tibble(
          peak_name = keep_peaks,
          group = label,
          host = cmat[keep_peaks, host_cols[idx]],
          ko   = cmat[keep_peaks, ko_cols[idx]]
        )
    }

    df <- dplyr::bind_rows(
        get_mean_df(1:5,   "Late ectoderm"),
        get_mean_df(6:8,   "Late mesoderm"),
        get_single_df(9,   "Early ectoderm"),
        get_mean_df(10:14, "Early mesoderm")
    ) %>%
    dplyr::mutate(
        group = factor(group, c("Early ectoderm", "Early mesoderm", "Late ectoderm", "Late mesoderm"))
    )

    stats_df <- df %>%
    dplyr::group_by(group) %>%
    dplyr::summarise(
      pct_gt1  = mean((ko - host) >  1, na.rm = TRUE) * 100,
      pct_lt_1 = mean((ko - host) < -1, na.rm = TRUE) * 100,
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      lab_tl = sprintf("%.1f%%", pct_gt1),
      lab_br = sprintf("%.1f%%", pct_lt_1),
      x_tl = -Inf, y_tl =  Inf,
      x_br =  Inf, y_br = -Inf
    )

    ggplot(df, aes(x = host, y = ko, color = group)) +
        geom_point(alpha = point_alpha, size = point_size) +
        geom_abline(slope = 1, intercept = 0, size = 0.4) +
        geom_abline(slope = 1, intercept = -1, linetype = "dashed", size = 0.4) +
        geom_abline(slope = 1, intercept = 1,  linetype = "dashed", size = 0.4) +
        facet_wrap(~group, nrow = 1) +
        scale_color_manual(values = c(
          "Late ectoderm"  = "#84a98c",
          "Late mesoderm"  = "#a594f9",
          "Early ectoderm" = "#b0c4b1",
          "Early mesoderm" = "#cdc1ff"
        )) +
        geom_text(
          data = stats_df,
          ggplot2::aes(x = x_tl, y = y_tl, label = lab_tl),
          inherit.aes = FALSE,
          hjust = -0.05, vjust = 1.3,
          size = 4.2, color = "black"
        ) +
        geom_text(
          data = stats_df,
          aes(x = x_br, y = y_br, label = lab_br),
          inherit.aes = FALSE,
          hjust = 1.05, vjust = -0.3,
          size = 4.2, color = "black"
        ) +
        theme_minimal(base_size = 20) +
        xlab("Host") +
        ylab("DKO") +
        theme(
          axis.text = ggplot2::element_text(size = 12),
          panel.border = ggplot2::element_rect(color = "black", fill = NA, size = 0.4),
          strip.text = ggplot2::element_text(face = "bold"),
          legend.position = "none"
        )
    }

plot_paired_box_by_cgbin <- function(
  cmat,
  peaks_df,
  name,
  clusters_keep,
  types_all,
  cluster_assignments_rds = "./data/output/%s_cluster_assignments.RDS",
  facet_nrow = 1,
  point_alpha = 0.15,
  line_alpha = 0.10,
  ymax = NULL,              # ignored; ylim fixed below
  y_bracket = -14.5,        # fixed bracket height
  bracket_drop = 0.15,      # how far the bracket "legs" drop
  y_star = -14.15,          # fixed star height
  p_adjust = "BH",          # "BH", "bonferroni", or "none"
  left_widths = c(1, 12)     # count panel vs box panel
) {

  # --- Choose peaks ---
  ann_rows <- readRDS(in_path(sprintf(cluster_assignments_rds, name)))
  cl_full <- suppressWarnings(as.integer(as.character(ann_rows$cluster)))
  if (all(is.na(cl_full))) cl_full <- as.integer(ann_rows$cluster)

  cl <- clusters_keep
  if (is.character(cl)) cl <- as.integer(gsub("^C", "", cl))
  cl <- as.integer(cl)

  keep_peaks <- intersect(rownames(ann_rows)[cl_full %in% cl], rownames(cmat))
  if (!length(keep_peaks)) stop("No matching peaks.")

  # --- cg_bin mapping ---
  peak_meta <- peaks_df %>%
    dplyr::select(peak_name, cg_bin) %>%
    dplyr::distinct() %>%
    dplyr::filter(peak_name %in% keep_peaks)

  if (!nrow(peak_meta)) stop("No cg_bin info for keep_peaks in peaks_df.")

  # --- columns ---
  host_cols <- vapply(types_all, `[[`, "", "host_col")
  ko_cols   <- vapply(types_all, `[[`, "", "ko_col")

  # --- same 4 group definitions ---
  get_mean_df <- function(idx, label) {
    tibble::tibble(
      peak_name = keep_peaks,
      group = label,
      host = rowMeans(cmat[keep_peaks, host_cols[idx], drop = FALSE], na.rm = TRUE),
      ko   = rowMeans(cmat[keep_peaks, ko_cols[idx],   drop = FALSE], na.rm = TRUE)
    )
  }
  get_single_df <- function(idx, label) {
    tibble::tibble(
      peak_name = keep_peaks,
      group = label,
      host = cmat[keep_peaks, host_cols[idx]],
      ko   = cmat[keep_peaks, ko_cols[idx]]
    )
  }

  wide_df <- dplyr::bind_rows(
    get_mean_df(1:5,   "Late ectoderm"),
    get_mean_df(6:8,   "Late mesoderm"),
    get_single_df(9,   "Early ectoderm"),
    get_mean_df(10:14, "Early mesoderm")
  ) %>%
    dplyr::inner_join(peak_meta, by = "peak_name") %>%
    dplyr::filter(!is.na(cg_bin), !is.na(host), !is.na(ko))

  group_levels <- c("Early ectoderm","Early mesoderm","Late ectoderm","Late mesoderm")
  wide_df$group <- factor(wide_df$group, levels = group_levels)

  # Long form for boxplots
  long_df <- wide_df %>%
    tidyr::pivot_longer(c(host, ko), names_to = "cond", values_to = "val") %>%
    dplyr::mutate(
      cond  = factor(cond, levels = c("host", "ko")),
      cg_bin = factor(cg_bin)
    )

  # --- Count bar (ALL peaks by cg_bin; matches function 1 style) ---
  count_df <- peak_meta %>%
    dplyr::group_by(cg_bin) %>%
    dplyr::summarise(n_peaks = dplyr::n_distinct(peak_name), .groups = "drop") %>%
    dplyr::mutate(
      prop = n_peaks / sum(n_peaks),
      x = "All peaks"
    ) %>%
    dplyr::mutate(cg_bin = factor(cg_bin, c("<0.02", "0.02-0.04", ">0.04")))

  p_count <- ggplot2::ggplot(count_df, ggplot2::aes(x = x, y = prop, fill = cg_bin)) +
    ggplot2::geom_col(width = 0.6) +
    ggplot2::scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
    ggplot2::scale_fill_manual(values = c("<0.02" = "grey70", "0.02-0.04" = "grey40", ">0.04" = "grey10")) +
    ggplot2::labs(
      x = NULL,
      y = "Fraction of peaks",
      fill = "CpG bin"
    ) +
    ggplot2::theme_minimal(base_size = 15) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(color = "black", fill = NA, size = 0.4),
      legend.position = "none"
    ) + theme(plot.margin = margin(5, 2, 5, 2)) + theme(
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank()
)

  # --- Paired Wilcoxon per (group, cg_bin) ---
  p_df <- wide_df %>%
    dplyr::group_by(group, cg_bin) %>%
    dplyr::summarise(
      n = sum(!is.na(host) & !is.na(ko)),
      p = dplyr::if_else(
        n >= 2,
        suppressWarnings(stats::wilcox.test(host, ko, paired = TRUE, exact = FALSE)$p.value),
        NA_real_
      ),
      .groups = "drop"
    )

  if (!identical(p_adjust, "none")) {
    p_df$p_adj <- stats::p.adjust(p_df$p, method = p_adjust)
  } else {
    p_df$p_adj <- p_df$p
  }

  p_df <- p_df %>%
    dplyr::mutate(
      star = dplyr::case_when(
        is.na(p_adj)   ~ "",
        p_adj <= 0.001 ~ "***",
        p_adj <= 0.01  ~ "**",
        p_adj <= 0.05  ~ "*",
        TRUE           ~ "ns"
      )
    )

  # --- bracket geometry coords (numeric x + offsets) ---
  x_levels <- levels(long_df$cg_bin)
  p_df <- p_df %>%
    dplyr::mutate(
      x_mid = as.numeric(factor(cg_bin, levels = x_levels)),
      x0    = x_mid - 0.22,
      x1    = x_mid + 0.22,
      y_br  = y_bracket,
      y_leg = y_bracket - bracket_drop,
      y_st  = y_star
    )

  # --- paired boxplot panel ---
  p_box <- ggplot2::ggplot(long_df, ggplot2::aes(x = cg_bin, y = val, fill = cond)) +
    ggplot2::geom_boxplot(
      position = ggplot2::position_dodge(width = 0.5),
      width = 0.4,
      outlier.shape = NA
    ) +
    ggplot2::facet_wrap(~group, nrow = facet_nrow) +
    ggplot2::labs(
      x = "CpG bin",
      y = "Accessibility",
      fill = NULL
    ) +
    ggplot2::scale_fill_manual(values = c("host" = "grey60", "ko" = "#a7c957"),
                               labels = c("host" = "Host", "ko" = "KO")) +
    ggplot2::theme_minimal(base_size = 15) +
    ggplot2::theme(
      axis.text.x  = ggplot2::element_text(size = 12),
      strip.text   = ggplot2::element_text(face = "bold"),
      panel.border = ggplot2::element_rect(color = "black", fill = NA, size = 0.4),
      legend.position = "top"
    ) +

    # ---- bracket ----
    ggplot2::geom_segment(
      data = p_df,
      ggplot2::aes(x = x0, xend = x0, y = y_br, yend = y_leg),
      inherit.aes = FALSE
    ) +
    ggplot2::geom_segment(
      data = p_df,
      ggplot2::aes(x = x1, xend = x1, y = y_br, yend = y_leg),
      inherit.aes = FALSE
    ) +
    ggplot2::geom_segment(
      data = p_df,
      ggplot2::aes(x = x0, xend = x1, y = y_br, yend = y_br),
      inherit.aes = FALSE
    ) +

    # ---- stars ----
    ggplot2::geom_text(
      data = p_df,
      ggplot2::aes(x = x_mid, y = y_st, label = star),
      inherit.aes = FALSE,
      size = 5
    ) +
    ggplot2::coord_cartesian(ylim = c(-20, -14))+ theme(
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank()
)

  combo <- (p_count | p_box) + patchwork::plot_layout(widths = left_widths)

  combo
}
