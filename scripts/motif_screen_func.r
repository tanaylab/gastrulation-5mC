library("parallel")

add_peak_energ = function(mod)
{

   npeaks = mod$npeaks
   npeaks_int = npeaks
   npeaks_int$start = (npeaks$start+npeaks$end)/2-150
   npeaks_int$end = (npeaks$start+npeaks$end)/2+150

   load("data/ref/motif_db.rda")
   motif_db_s <- motif_db %>% filter(dataset %in% c("HOMER", "JASPAR", "JOLMA") | motif %in% c("SCENIC.cisbp__M00615", "SCENIC.cisbp__M00626"))

   if(file.exists("data/output/npeaks_energ_raw.RDS")) {
      npeaks_energ = readRDS(in_path("data/output/npeaks_energ_raw.RDS"))
   } else {
      message("will compute energies")
#normalization in IQ is problematic
#     npeaks_energ <- compute_motif_energies(peak_intervals=npeaks_int,
#                 motif_db_s, normalization_intervals = norm_intervs)

      npeaks_energ <- prego::gextract_pwm(npeaks_int, dataset = motif_db_s, prior = 0.01)

      saveRDS(npeaks_energ, out_path("data/output/npeaks_energ_raw.RDS"))
   }
   mod$npeaks_energ = npeaks_energ
   return(mod)
}

stratify_stat_by_pwms = function(mod, energ, bin, st_abs, st_diff, 
											peak_filter, back_percentile=0.98)
{
	maxbin = max(bin, na.rm=T)
	test_e_linkage = function(i, f, st_abs, st_diff)
	{
		x = energ[f,i]
		y = st_diff[f]
		y_abs = st_abs[f,1]
	#	et=quantile(x,0.95, na.rm=T); 
	#	em=quantile(x,0.5, na.rm=T); 
		et=quantile(energ[peak_filter, i], back_percentile, na.rm=T); 
		em=quantile(energ[peak_filter, i], 0.95, na.rm=T); 
		n_hit = sum(x >= et, na.rm=T)
		n_low = sum(x <= em, na.rm=T)
		n = sum(f,na.rm=T)
		if(n_hit > 100 & n_low > 100) {
		#thus in the two-sample case 'alternative = greater includes distributions for which 'x' is stochastically _smaller_ than 'y'
			s_abs = ks.test(y_abs[x>=et],y_abs[x<=em], alternative="less")
			s_up = ks.test(y[x>=et],y[x<=em], alternative="less") 
			s_down = ks.test(y[x>=et],y[x<=em], alternative="greater") 
			return(c(s_up$statistic, s_up$p.value, s_down$statistic, s_down$p.value, s_abs$statistic, s_abs$p.value, n_hit, n))
		} else {
			return(rep(NA, 8))
		}
	}

	strat_highmot = function(i, st_abs, min_n = 10) {
		x = energ[peak_filter,i]
		et=quantile(energ[peak_filter, i], back_percentile, na.rm=T); 
		y_abs = st_abs[peak_filter,]
		trend1 = tapply(y_abs[x >= et, 1], bin[peak_filter][x >= et], mean, na.rm=T)
		trend1_med = tapply(y_abs[x >= et, 1], bin[peak_filter][x >= et], median,na.rm=T)
		trend2 = tapply(y_abs[x >= et, 2], bin[peak_filter][x >= et], mean,na.rm=T)
		trend2_med = tapply(y_abs[x >= et, 2], bin[peak_filter][x >= et], median, na.rm=T)
		trend4 = tapply(y_abs[x >= et, 4], bin[peak_filter][x >= et], mean, na.rm=T)
		trend4_med = tapply(y_abs[x >= et, 4], bin[peak_filter][x >= et], median, na.rm=T)
		trend_n = tapply(y_abs[x >= et, 1], bin[peak_filter][x >= et], length)
		trend1[trend_n < min_n] = NA
		trend1_med[trend_n < min_n] = NA
		trend2[trend_n < min_n] = NA
		trend2_med[trend_n < min_n] = NA
		trend4[trend_n < min_n] = NA
		trend4_med[trend_n < min_n] = NA
		allbins = as.character(1:maxbin)
		return(c(trend1[allbins], trend2[allbins], trend4[allbins],
					trend1_med[allbins], trend2_med[allbins], trend4_med[allbins],
					trend_n[allbins]))
	}

	trends = mclapply(1:ncol(energ), strat_highmot, st_abs = st_abs, mc.cores=40)
	trends = do.call('rbind', trends)

	ds_up = NULL
	ds_down = NULL
	ds_abs = NULL
	ds_up_p = NULL
	ds_down_p = NULL
	ds_abs_p = NULL
	ds_hit = NULL

	for(i in 1:maxbin) {
		message("i is ",i)
		f = (bin==i) & peak_filter

		ds = mclapply(1:ncol(energ), test_e_linkage, f, st_abs = st_abs, st_diff = st_diff, mc.cores=40)
		ds_df = do.call('cbind', ds)
	
		if(is.null(ds_up)) {	
			ds_up = t(as.matrix(ds_df[1,]))
			ds_up_p = t(as.matrix(ds_df[2,]))
			ds_down = t(as.matrix(ds_df[3,]))
			ds_down_p = t(as.matrix(ds_df[4,]))
			ds_abs = t(as.matrix(ds_df[5,]))
			ds_abs_p = t(as.matrix(ds_df[6,]))
			ds_hit = t(as.matrix(ds_df[7,]))
		} else {
			ds_up = rbind(ds_up, ds_df[1,])
			ds_up_p = rbind(ds_up_p, ds_df[2,])
			ds_down = rbind(ds_down, ds_df[3,])
			ds_down_p = rbind(ds_down_p, ds_df[4,])
			ds_abs = rbind(ds_abs, ds_df[5,])
			ds_abs_p = rbind(ds_abs_p, ds_df[6,])
			ds_hit = rbind(ds_hit, ds_df[7,])
		}
	}
	colnames(ds_up) = colnames(energ)
	colnames(ds_up_p) = colnames(energ)
	colnames(ds_down) = colnames(energ)
	colnames(ds_down_p) = colnames(energ)
	colnames(ds_abs) = colnames(energ)
	colnames(ds_abs_p) = colnames(energ)
	colnames(ds_hit) = colnames(energ)
	rownames(trends) = colnames(energ)
	return(list(trends = trends, ds_up = ds_up, ds_up_p = ds_up_p, 
					ds_down = ds_down, ds_down_p = ds_down_p,
					ds_abs = ds_abs, ds_abs_p = ds_abs_p, ds_hit = ds_hit))
}

add_pwms_stats_TSS = function(mod, abs_T = 0)
{
	if(file.exists("data/output/pwms_stats_TSS.RDS")) {
		mod$pwms_stats_TSS = readRDS(in_path("data/output/pwms_stats_TSS.RDS"))
		return(mod)
	}
	energ = mod$npeaks_energ[,-c(1:3,6:8)]

	cgbin = as.numeric(cut(mod$npeaks$CG, c(0,seq(0.02,0.08,by = 0.02),1), include.lowest = T))

	lp = mod$lin_peaks

	u_acp = lp$acp_ds
	e_acp  = t(t(u_acp)/colSums(u_acp))
	le_acp = log2(1e-6+e_acp)
	acp_lf = le_acp[,2]-le_acp[,4]

	peak_filter = abs(mod$npeaks$d_tss)<=1000 & le_acp[,1] < abs_T
	pwms_cg_acp = stratify_stat_by_pwms(mod, energ = energ, bin=cgbin, st_abs=le_acp, st_diff=acp_lf, peak_filter = peak_filter)

	u_rnp = lp$rnp_ds
	e_rnp  = t(t(u_rnp)/colSums(u_rnp))
	le_rnp = log2(1e-6+e_rnp)
	rnp_lf = le_rnp[,2]-le_rnp[,4]

	peak_filter = abs(mod$npeaks$d_tss)<=1000 & le_rnp[,1] < abs_T
	pwms_cg_rnp = stratify_stat_by_pwms(mod, energ = energ, bin=cgbin, st_abs=le_rnp, st_diff=rnp_lf, peak_filter = peak_filter)

	u_meso = lp$meso_ds
	e_meso  = t(t(u_meso)/colSums(u_meso))
	le_meso = log2(1e-6+e_meso)
	meso_lf = le_meso[,2]-le_meso[,4]

	peak_filter = abs(mod$npeaks$d_tss)<=1000 & le_meso[,1] < abs_T
	pwms_cg_meso = stratify_stat_by_pwms(mod, energ = energ, bin=cgbin, st_abs=le_meso, st_diff=meso_lf, peak_filter = peak_filter)

	pwms_stats_TSS = list(acp = pwms_cg_acp, rnp = pwms_cg_rnp, meso = pwms_cg_meso)
	saveRDS(pwms_stats_TSS, out_path("data/output/pwms_stats_TSS.RDS"))
	mod$pwms_stats_TSS = pwms_stats_TSS

	return(mod)
}

add_pwms_stats_v_wt = function(mod, abs_T = 0)
{
	if(file.exists("data/output/pwms_stat_v_wt.RDS")) {
		mod$pwms_stat_v_wt = readRDS(in_path("data/output/pwms_stat_v_wt.RDS"))
		return(mod)
	}
	energ = mod$npeaks_energ[,-c(1:3,6:8)]

	cgbin = as.numeric(cut(mod$npeaks$CG, c(seq(0,0.05,0.01),1)))

	lp = mod$lin_peaks

	u_acp = lp$acp_ds
	e_acp  = t(t(u_acp)/colSums(u_acp))
	le_acp = log2(1e-6+e_acp)

	# DNMT - WT
	acp_lf = le_acp[,2]-le_acp[,1]
	peak_filter = abs(mod$npeaks$d_tss)>1000 & le_acp[,1] < abs_T
	dnmt_pwms_cg_acp = stratify_stat_by_pwms(mod, energ = energ, bin=cgbin, st_abs=le_acp, st_diff=acp_lf, peak_filter = peak_filter)
	# WT - TET
	acp_lf = le_acp[,1]-le_acp[,4]
	peak_filter = abs(mod$npeaks$d_tss)>1000 & le_acp[,1] < abs_T
	tet_pwms_cg_acp = stratify_stat_by_pwms(mod, energ = energ, bin=cgbin, st_abs=le_acp, st_diff=acp_lf, peak_filter = peak_filter)

	u_rnp = lp$rnp_ds
	e_rnp  = t(t(u_rnp)/colSums(u_rnp))
	le_rnp = log2(1e-6+e_rnp)

	# DNMT - WT
	rnp_lf = le_rnp[,2]-le_rnp[,1]
	peak_filter = abs(mod$npeaks$d_tss)>1000 & le_rnp[,1] < abs_T
	dnmt_pwms_cg_rnp = stratify_stat_by_pwms(mod, energ = energ, bin=cgbin, st_abs=le_rnp, st_diff=rnp_lf, peak_filter = peak_filter)
	# WT - TET
	rnp_lf = le_rnp[,1]-le_rnp[,4]
	peak_filter = abs(mod$npeaks$d_tss)>1000 & le_rnp[,1] < abs_T
	tet_pwms_cg_rnp = stratify_stat_by_pwms(mod, energ = energ, bin=cgbin, st_abs=le_rnp, st_diff=rnp_lf, peak_filter = peak_filter)

	u_meso = lp$meso_ds
	e_meso  = t(t(u_meso)/colSums(u_meso))
	le_meso = log2(1e-6+e_meso)

	# DNMT - WT
	meso_lf = le_meso[,2]-le_meso[,1]
	peak_filter = abs(mod$npeaks$d_tss)>1000 & le_meso[,1] < abs_T
	dnmt_pwms_cg_meso = stratify_stat_by_pwms(mod, energ = energ, bin=cgbin, st_abs=le_meso, st_diff=meso_lf, peak_filter = peak_filter)
	# WT - TET
	meso_lf = le_meso[,1]-le_meso[,4]
	peak_filter = abs(mod$npeaks$d_tss)>1000 & le_meso[,1] < abs_T
	tet_pwms_cg_meso = stratify_stat_by_pwms(mod, energ = energ, bin=cgbin, st_abs=le_meso, st_diff=meso_lf, peak_filter = peak_filter)

	pwms_stat = list(dnmt_acp = dnmt_pwms_cg_acp, dnmt_rnp = dnmt_pwms_cg_rnp, dnmt_meso = dnmt_pwms_cg_meso,
	                 tet_acp = tet_pwms_cg_acp, tet_rnp = tet_pwms_cg_rnp, tet_meso = tet_pwms_cg_meso)
	saveRDS(pwms_stat, out_path("data/output/pwms_stat_v_wt.RDS"))
	mod$pwms_stat_v_wt = pwms_stat

	return(mod)
}

screen_pwms_TSS <- function(ps, T_ks = 0.25){
                        
    foc_pwms_rnp = find_sig_pwms_mod(ps$rnp)
    foc_pwms_acp = find_sig_pwms_mod(ps$acp)
    foc_pwms_meso = find_sig_pwms_mod(ps$meso)

    foc_pwms <- list(foc_pwms_rnp, foc_pwms_acp, foc_pwms_meso)
    ks_up = get_ks_mat(foc_pwms, dir = "u")
    up = names(which(apply(ks_up$maxks, 1, function(x) {max(x) > T_ks})))

    ks_down = get_ks_mat(foc_pwms, dir = "d")
    down = names(which(apply(ks_down$maxks, 1, function(x) {max(x) > T_ks})))

    max_ks = do.call(cbind, list(
        apply(ks_up$maxks, 1, max), 
        apply(ks_down$maxks, 1, max)))

    nms = union(up, down)
    ks = max_ks[nms,]
    colnames(ks) = c("u", "d")

    return(list(foc_pwms = nms, max_ks = ks, qvals = cbind(ks_up$qvals, ks_down$qvals)))
    
}

screen_pwms <- function(ps, T_ks = 0.25){
                        
    foc_pwms_rnp = find_sig_pwms_mod(ps$dnmt_rnp)
    foc_pwms_acp = find_sig_pwms_mod(ps$dnmt_acp)
    foc_pwms_meso = find_sig_pwms_mod(ps$dnmt_meso)

    dnmt_pwms <- list(foc_pwms_rnp, foc_pwms_acp, foc_pwms_meso)
    ks_dnmt_up = get_ks_mat(dnmt_pwms, dir = "u")
    dnmt_up = names(which(apply(ks_dnmt_up$maxks, 1, function(x) {max(x) > T_ks})))

    ks_dnmt_down = get_ks_mat(dnmt_pwms, dir = "d")
    dnmt_down = names(which(apply(ks_dnmt_down$maxks, 1, function(x) {max(x) > T_ks})))

    foc_pwms_rnp = find_sig_pwms_mod(ps$tet_rnp)
    foc_pwms_acp = find_sig_pwms_mod(ps$tet_acp)
    foc_pwms_meso = find_sig_pwms_mod(ps$tet_meso)

    tet_pwms <- list(foc_pwms_rnp, foc_pwms_acp, foc_pwms_meso)
    ks_tet_up = get_ks_mat(tet_pwms, dir = "u")
    tet_up = names(which(apply(ks_tet_up$maxks, 1, function(x) {max(x) > T_ks})))

    ks_tet_down = get_ks_mat(tet_pwms, dir = "d")
    tet_down = names(which(apply(ks_tet_down$maxks, 1, function(x) {max(x) > T_ks})))

    max_ks = do.call(cbind, list(
        apply(ks_dnmt_up$maxks, 1, max), 
        apply(ks_dnmt_down$maxks, 1, max), 
        apply(ks_tet_up$maxks, 1, max), 
        apply(ks_tet_down$maxks, 1, max)))

    nms = union(union(dnmt_up, dnmt_down), union(tet_up, tet_down))
    ks = max_ks[nms,]
    colnames(ks) = c("du", "dd", "tu", "td")

	all_qvals = cbind(ks_dnmt_up$qvals, ks_dnmt_down$qvals, ks_tet_up$qvals, ks_tet_down$qvals)
	rownames(all_qvals) = rownames(max_ks)
	colnames(all_qvals) = c(paste0("du_", 1:3), paste0("dd_", 1:3), paste0("tu_", 1:3), paste0("td_", 1:3))

    return(list(foc_pwms = nms, max_ks = ks, qvals = all_qvals))
    
}

   

# }

plot_clust_pwms <- function(ps, screen_res, suffix = NULL){
    
    nms = screen_res$foc_pwms
    ks = screen_res$ks
    
    tr_rnp = ps$dnmt_rnp$trends
    tr_acp = ps$dnmt_acp$trends
    tr_meso = ps$dnmt_meso$trends

    abs = cbind(tr_rnp[nms,1:6], tr_acp[nms,1:6], tr_meso[nms,1:6])
    diff = cbind(tr_acp[nms,1:6]-tr_rnp[nms,1:6], tr_meso[nms,1:6]-tr_rnp[nms,1:6])
    dnmt = cbind(tr_rnp[nms, 6+1:6]-tr_rnp[nms, 0+1:6], tr_acp[nms, 6+1:6]-tr_acp[nms, 0+1:6],tr_meso[nms, 6+1:6]-tr_meso[nms, 0+1:6])
    tet = cbind(tr_rnp[nms,0+1:6]-tr_rnp[nms,12+1:6], tr_acp[nms,0+1:6]-tr_acp[nms,12+1:6],tr_meso[nms,0+1:6]-tr_meso[nms,12+1:6])

	# Save mats to RDSs
	saveRDS(abs, out_path("data/output/pwm_trends_abs.RDS"))
	saveRDS(diff, out_path("data/output/pwm_trends_diff.RDS"))
	saveRDS(dnmt, out_path("data/output/pwm_trends_dnmt.RDS"))
	saveRDS(tet, out_path("data/output/pwm_trends_tet.RDS"))

    meth = calc_meth_stats(mod, nms)
    meth = meth[nms,]

    abs_n = t(t(abs)-colMeans(abs,na.rm=T))
    feats = cbind(abs_n, dnmt*2, tet*6)
    colnames(feats) = 1:ncol(feats)

    k = ifelse(nrow(feats) > 16, 16, 1)
    km = tglkmeans::TGL_kmeans(as.matrix(feats), k, seed=19)
    saveRDS(km, out_path("data/output/pwm_trends_km.RDS"))
    gaps_row = which(!duplicated(sort(km$cluster)))
    gaps_row = gaps_row[2:length(gaps_row)]
    gaps_row = gaps_row - 1
    if (k==1) {gaps_row = NULL}
    
    rshades = colorRampPalette(c("white","pink", "darkred", "purple"))(1000)
    bshades = colorRampPalette(c("white", "lightblue1", "lightblue2", "deepskyblue2", "deepskyblue4", "black"))(1000)
    gshades = colorRampPalette(c("white", "grey95", "darkolivegreen3", "darkolivegreen4", "darkolivegreen", "black"))(1000)
    mshades = colorRampPalette(c("white", "wheat", "goldenrod1", "darkorange", "saddlebrown", "black"))(1000)
    bwshades = colorRampPalette(c("white", "plum"))(1000)
        
    # DNMT - WT
    show_fig("figs/Fig5/final_heatmaps/pwm_km_dnmt.png", draw = pheatmap::pheatmap(dnmt[order(km$cluster),], cluster_rows=F, cluster_cols=F, col=bshades, 
                       w=2,h=15, legend = F, show_colnames=T, angle_col = 0, show_rownames=F, breaks=c(-0.2,seq(0.1,0.7,l=999),0.8), gaps_row = gaps_row, gaps_col = c(6,12)))
    # WT - TET
    show_fig("figs/Fig5/final_heatmaps/pwm_km_tet.png", draw = pheatmap::pheatmap(tet[order(km$cluster),], cluster_rows=F, cluster_cols=F, col=rshades, 
                       w=2,h=15, legend = F, show_colnames=T, angle_col = 0, show_rownames=F, breaks=c(-0.2,seq(0,0.35,l=999),0.5), gaps_row = gaps_row, gaps_col = c(6,12)))
    show_fig("figs/Fig5/final_heatmaps/pwm_km_trends.png", draw = pheatmap::pheatmap(abs[order(km$cluster),], cluster_rows=F, cluster_cols=F, col=gshades, 
                       w=2,h=15, legend = F, show_colnames=T, angle_col = 0, show_rownames=F, breaks=seq(-20,-16.5,l=999), gaps_row = gaps_row, gaps_col = c(6,12)))
    show_fig("figs/Fig5/final_heatmaps/pwm_km_meth.png", draw = pheatmap::pheatmap(meth[order(km$cluster),], cluster_rows=F, cluster_cols=F, col=mshades, 
                       w=3,h=15, legend = F, fontsize_row = 5, show_colnames=T, angle_col = 0, show_rownames=T, breaks=seq(0,1,l=999), gaps_row = gaps_row, gaps_col = c(6)))	
    
    
    # Plot scales
    marg = c(2,1,2,15)
    vec <- c(-0.2,seq(0.1,0.7,l=999),0.8)
    points <- seq(from = 1, to = length(vec), length.out = 7)
    selected_points <- vec[round(points)]
    show_fig("./figs/Fig5/scales_to_scale/5A_scales/top_dnmt_scale.pdf", show = FALSE, width = 7, height = 20, draw = {
        par(mar = marg)
        image(1, 1:100, t(as.matrix(1:100)), col = bshades, axes = FALSE, xlab = "", ylab = "")
        axis(4, at = seq(1, 100, length.out = 7), labels = round(selected_points, 2), las = 1, cex.axis = 5)
    })

    vec <- c(-0.2, seq(0, 0.35, l = 999), 0.5)
    points <- seq(from = 1, to = length(vec), length.out = 7)
    selected_points <- vec[round(points)]
    show_fig("./figs/Fig5/scales_to_scale/5A_scales/top_tet_scale.pdf", show = FALSE, width = 7, height = 20, draw = {
        par(mar = marg)
        image(1, 1:100, t(as.matrix(1:100)), col = rshades, axes = FALSE, xlab = "", ylab = "")
        axis(4, at = seq(1, 100, length.out = 7), labels = round(selected_points, 2), las = 1, cex.axis = 5)
    })

    vec <- seq(-20,-16.5,l=999)
    points <- seq(from = 1, to = length(vec), length.out = 7)
    selected_points <- vec[round(points)]
    show_fig("./figs/Fig5/scales_to_scale/5A_scales/top_wt_scale.pdf", show = FALSE, width = 7, height = 20, draw = {
        par(mar = marg)
        image(1, 1:100, t(as.matrix(1:100)), col = gshades, axes = FALSE, xlab = "", ylab = "")
        axis(4, at = seq(1, 100, length.out = 7), labels = round(selected_points, 2), las = 1, cex.axis = 5)
    })

    vec <- seq(0,1,l=999)
    points <- seq(from = 1, to = length(vec), length.out = 7)
    selected_points <- vec[round(points)]
    show_fig("./figs/Fig5/scales_to_scale/5A_scales/top_meth_scale.pdf", show = FALSE, width = 7, height = 20, draw = {
        par(mar = marg)
        image(1, 1:100, t(as.matrix(1:100)), col = mshades, axes = FALSE, xlab = "", ylab = "")
        axis(4, at = seq(1, 100, length.out = 7), labels = round(selected_points, 2), las = 1, cex.axis = 5)
    })
    
}

plot_clust_pwms_TSS <- function(ps, screen_res){
    
    nms = screen_res$foc_pwms
    ks = screen_res$ks
    
    tr_rnp = ps$rnp$trends
    tr_acp = ps$acp$trends
    tr_meso = ps$meso$trends

    abs = cbind(tr_rnp[nms,1:5], tr_acp[nms,1:5], tr_meso[nms,1:5])
	colnames(abs) = rep(c(1:5), 3)
    diff = cbind(tr_acp[nms,1:5]-tr_rnp[nms,1:5], tr_meso[nms,1:5]-tr_rnp[nms,1:5])
    dnmt = cbind(tr_rnp[nms, 5+1:5]-tr_rnp[nms, 0+1:5], tr_acp[nms, 5+1:5]-tr_acp[nms, 0+1:5],tr_meso[nms, 5+1:5]-tr_meso[nms, 0+1:5])
	colnames(dnmt) = rep(c(1:5), 3)
    tet = cbind(tr_rnp[nms,0+1:5]-tr_rnp[nms,10+1:5], tr_acp[nms,0+1:5]-tr_acp[nms,10+1:5],tr_meso[nms,0+1:5]-tr_meso[nms,10+1:5])
	colnames(tet) = rep(c(1:5), 3)
    meth = calc_meth_stats(mod, nms, cgrange = c(0,seq(0.02,0.08,by = 0.02),1))
    meth = meth[nms,]

    abs_n = t(t(abs)-colMeans(abs,na.rm=T))
    feats = cbind(abs_n, dnmt*2, tet*6)
    colnames(feats) = 1:ncol(feats)

    k = ifelse(nrow(feats) > 16, 16, 1)
    km = tglkmeans::TGL_kmeans(as.matrix(feats), k, seed=19)
    saveRDS(km, out_path("data/output/pwm_trends_km.RDS"))
    gaps_row = which(!duplicated(sort(km$cluster)))
    gaps_row = gaps_row[2:length(gaps_row)]
    gaps_row = gaps_row - 1
    if (k==1) {gaps_row = NULL}
    
    rshades = colorRampPalette(c("white","pink", "darkred", "purple"))(1000)
    bshades = colorRampPalette(c("white", "lightblue1", "lightblue2", "deepskyblue2", "deepskyblue4", "black"))(1000)
    gshades = colorRampPalette(c("white", "grey95", "darkolivegreen3", "darkolivegreen4", "darkolivegreen", "black"))(1000)
    mshades = colorRampPalette(c("white", "wheat", "goldenrod1", "darkorange", "saddlebrown", "black"))(1000)
    bwshades = colorRampPalette(c("white", "plum"))(1000)
        
    # DNMT - WT
    show_fig("figs/Fig5/final_heatmaps/pwm_km_dnmt_TSS.pdf", draw = pheatmap::pheatmap(dnmt[order(km$cluster),], cluster_rows=F, cluster_cols=F, col=bshades, 
                       w=2,h=5, legend = F, show_colnames=T, angle_col = 0, show_rownames=F, breaks=c(-0.2,seq(0.1,0.7,l=999),0.8), gaps_row = gaps_row, gaps_col = c(5,10)))
    # WT - TET
    show_fig("figs/Fig5/final_heatmaps/pwm_km_tet_TSS.pdf", draw = pheatmap::pheatmap(tet[order(km$cluster),], cluster_rows=F, cluster_cols=F, col=rshades, 
                       w=2,h=5, legend = F, show_colnames=T, angle_col = 0, show_rownames=F, breaks=c(-0.2,seq(0,0.35,l=999),0.5), gaps_row = gaps_row, gaps_col = c(5,10)))
    show_fig("figs/Fig5/final_heatmaps/pwm_km_trends_TSS.pdf", draw = pheatmap::pheatmap(abs[order(km$cluster),], cluster_rows=F, cluster_cols=F, col=gshades, 
                       w=2,h=5, legend = F, show_colnames=T, angle_col = 0, show_rownames=F, breaks=seq(-20,-16.5,l=999), gaps_row = gaps_row, gaps_col = c(5,10)))
    show_fig("figs/Fig5/final_heatmaps/pwm_km_meth_TSS.pdf", draw = pheatmap::pheatmap(meth[order(km$cluster),], cluster_rows=F, cluster_cols=F, col=mshades, 
                       w=3,h=5, legend = F, fontsize_row = 5, show_colnames=T, angle_col = 0, show_rownames=T, breaks=seq(0,1,l=999), gaps_row = gaps_row, gaps_col = c(5)))	
}

find_sig_pwms_mod = function(pstat, T_ks = 0.25, T_pv = 1e-4)
{
	dsu = pstat$ds_up[1:5,]
	pvu = pstat$ds_up_p[1:5,]

	dsd = pstat$ds_down[1:5,]
	pvd = pstat$ds_down_p[1:5,]

	dsu[is.na(dsu)] = 0
	pvu[is.na(pvu)] = 1

	dsd[is.na(dsd)] = 0
	pvd[is.na(pvd)] = 1
    
    return(list("dsu" = dsu, "pvu" = pvu, "dsd" = dsd, "pvd" = pvd))
}

extract_stats <- function(pwms, dir){
    
    if (dir == "u") {
        uvals <- pwms$dsu
        pvals <- pwms$pvu
    } else if (dir == "d") {
        uvals <- pwms$dsd
        pvals <- pwms$pvd
    }
    
    pvals_maxks = rep(0, ncol(uvals))
    for (i in 1:ncol(uvals)){
        j = which.max(uvals[,i])
        pvals_maxks[i] <- pvals[j,i]                    
    }
    
    maxks = apply(uvals, 2, max)
    
    return(list("maxks" = maxks, "pvals_maxks" = pvals_maxks))
}

get_ks_mat <- function(pwms, dir, T_pv = 1e-4){
    
    results <- lapply(pwms, extract_stats, dir = dir)

    maxks_combined <- do.call(cbind, lapply(results, function(res) res$maxks))
    pvals_maxks_combined <- do.call(cbind, lapply(results, function(res) res$pvals_maxks))

    maxks_combined = maxks_combined[3:nrow(maxks_combined),] 
    pvals_maxks_combined = pvals_maxks_combined[3:nrow(pvals_maxks_combined),]
	message("nrow maxks_combined: ", nrow(maxks_combined))
	message("ncol maxks_combined: ", ncol(maxks_combined))

	# FDR correction
	pvals_vector <- as.vector(pvals_maxks_combined)
	pvals_fdr_vector <- p.adjust(pvals_vector, method = "fdr")
	qvals_matrix <- matrix(pvals_fdr_vector, nrow = nrow(pvals_maxks_combined), ncol = ncol(pvals_maxks_combined))                                        
	rownames(qvals_matrix) <- rownames(pvals_maxks_combined)

    maxks_combined[qvals_matrix > T_pv] = 0
                                   
    return(list("maxks" = maxks_combined, "qvals" = qvals_matrix))
}

calc_meth_stats = function(mod, foc_pwms, cgrange = c(seq(0,0.05,0.01),1)){
	energ = mod$npeaks_energ[,-c(1:3,6:8)]
	npeaks = mod$npeaks
	npeaks$cgbin = as.numeric(cut(npeaks$CG, cgrange))
	peak_filter = abs(mod$npeaks$d_tss)>1000

	e_foc = energ[peak_filter, foc_pwms]
	tope_mask = apply(e_foc, 2, function(x) {x > quantile(x, 0.98, na.rm=T)})

	npeaks_f = npeaks[peak_filter,]
	m = gextract_meth(c("CHIM.comb_wt_E8_5", "CHIM.comb_tet_tko_E8_5"), names = c("wt_m", "tet_m"), 
						intervals = npeaks_f, iterator = npeaks_f, min_cov = 10)

	npeaks_f$wt_m = m$wt_m
	npeaks_f$tet_m = m$tet_m

	cg_m_plus_masks = cbind(npeaks_f, tope_mask) %>% filter(!is.na(cgbin))

	mean_by_cg_bin <- function(df, selected_cols, stat_col) {
		map_dfr(selected_cols, ~ {
			df %>%
			filter(.data[[.x]]) %>%
			group_by(cgbin) %>%
			summarize(mean = mean(.data[[stat_col]], na.rm = TRUE), .groups = 'drop') %>%
			mutate(motif = .x)
			}
		)
	}

	long_wt = mean_by_cg_bin(cg_m_plus_masks, foc_pwms, "wt_m")
	wt_m = long_wt %>% pivot_wider(id_cols = "motif", names_from = "cgbin", values_from = "mean") %>%
		column_to_rownames("motif")

	long_tet = mean_by_cg_bin(cg_m_plus_masks, foc_pwms, "tet_m")
	tet_m = long_tet %>% pivot_wider(id_cols = "motif", names_from = "cgbin", values_from = "mean") %>%
		column_to_rownames("motif")

	motif_meth_trends = cbind(wt_m, tet_m)
	return(motif_meth_trends)
}
