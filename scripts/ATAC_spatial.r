library("data.table")

get_cmp_mats = function(mod, state, intervs, dnmt_dko_mod = "pool", window=1000, smooth_w = 0)
{
	if(!dnmt_dko_mod %in% c("pool", "ignore", "compare")) {
		message("dko mod should be one of pool, ignore, compare")
		return(NULL)
	}

	tr_exps =  c()
	tr_karyo = list()
	tr_dsamp_r = c()
	gen_track_vext(state, w_ext = 0)
	tr_exps["wt"] = "rowSums(cbind(wtl_ext, wts_ext, wtht_ext, wthd_ext, wthtet_ext),na.rm=T)"
	tr_karyo[["wt"]] = NULL

	if(dnmt_dko_mod == "ignore") {
		tr_exps["dnmt"] = "dnmt_ext"
		tr_karyo[["dnmt"]] = mod$r_chr_cov_dnmt[state,]
	} else if(dnmt_dko_mod == "pool") {
		tr_exps["dnmt"] = "rowSums(cbind(dnmt_ext, dnmtd_ext),na.rm=T)"
		tr_karyo[["dnmt"]] = mod$r_chr_cov_dnmt_both[state,]
	} else {
		tr_exps["dnmt"] = "dnmt_ext"
		tr_karyo[["dnmt"]] = mod$r_chr_cov_dnmt[state,]
		tr_exps["dnmtd"] = "dnmtd_ext"
		tr_karyo[["dnmtd"]] = mod$r_chr_cov_dnmtd[state,]
	}
	tr_exps["tet"] = "tet_ext"
	tr_karyo[["tet"]] = mod$r_chr_cov_tet[state,]

	tr_dsamp_r["wt"] = mod$r_to_dsamp[[dnmt_dko_mod]][state, "wt"]
	tr_dsamp_r["dnmt"] = mod$r_to_dsamp[[dnmt_dko_mod]][state, "dnmt"]
	if(dnmt_dko_mod == "compare") {
		tr_dsamp_r["dnmtd"] = mod$r_to_dsamp[[dnmt_dko_mod]][state, "dnmtd"]
	}
	tr_dsamp_r["tet"] = mod$r_to_dsamp[[dnmt_dko_mod]][state, "tet"]

	spat = get_spat_mat(tr_exps = tr_exps, 	
				tr_karyo = tr_karyo, 
				tr_dsamp_r = tr_dsamp_r, 
 				intervs = intervs, w = window,
				smooth_w = smooth_w)
	
	return(spat)	
}

get_spat_mat = function(tr_exps, tr_karyo, tr_dsamp_r, intervs, w=1000, smooth_w=0)
{
	tr_nms = names(tr_exps)
	center = intervs

   center$start = center$center
   center$end = center$center + 1
    
	gvtrack.create("d_center", center, "distance")

	region = center[,c("chrom","start","end","strand")]
	region$start = center$start - w
	region$end = center$end + (w + 20)
		
	spat = gextract(c("d_center", tr_exps),
			intervals = region, 
			iterator=20, colnames=c("d_center", tr_nms))

	spat[is.na(spat)] = 0
	message("into karyo correct")
	for(nm in tr_nms) {
		if(!is.null(tr_karyo[[nm]])) {
			spat = downsamp_karyo_balance(spat, nm, tr_karyo[[nm]])
		}
	}
	for(nm in tr_nms) {
		if(!is.null(tr_dsamp_r[nm])) {
			to_ds = matrix(spat[,nm],ncol=1)
			r_ds = tr_dsamp_r[nm]
			message("into dsamp, ", nm, " targ ", floor(sum(to_ds)*r_ds))
			post_ds = downsample_matrix(to_ds,  target_n = floor(sum(to_ds)*r_ds), seed=19)
			spat[,nm] = post_ds
		}
	}
#karyo dsamp
#cov dsampe

	spat$d_center_bin = floor(spat$d_center/20)
	maxbin = w/20

	mats = list()
	for(nm in tr_nms) {
		message("generate wide format for ", nm)
		mat = spat[abs(spat$d_center_bin)<=maxbin,] %>% reshape2::dcast(intervalID ~ d_center_bin, value.var=nm, fun.aggregate=mean,fill=0)
		mat[is.na(mat)] = 0
		mat  = as.matrix(mat[,-1])
		if(smooth_w > 1) {
			mat = t(apply(mat, 1, function(x) rollmean(x,smooth_w,f='e')))
		}
		mats[[nm]] = mat
	}

	return(mats)
}

reord_smooth_on_km = function(mat, km, w_size, should_hc = T)
{
    for(i in 1:max(km$cluster)){
		clmat = mat[km$cluster == i,]
		if(should_hc) {
			hc = hclust(tgs_dist(log2(0.25+clmat)), "ward.D2")
			clmat = clmat[hc$order,]
		}
		clmat_sm = apply(clmat, 2, function(x) rollmean(x, w_size, na.rm=T, f='e'))
		if(i==1) {
			mat_sm = clmat_sm
		} else {
			mat_sm = rbind(mat_sm, clmat_sm)
		}
	}
	return(mat_sm)
}

# for each we cluster atac on wt/dnmt/tet X conditions. Showing heat map, and also profile per clust (compare to avg)
cluster_ehs_spat_mats = function(mod, modes, states, max_k = 32, min_n_per_km=200, km_seed = NULL, filter_low_signal = F)
{
	wn = 101

	c_nms = c("WT_lg.tot_rp", "WT_sm.tot_rp", "DNMT_TKO_host.tot_rp", "DNMT_DKO_host.tot_rp", "TET_TKO_host.tot_rp")
	tot_wt_rp = rowSums(mod$cov_stat[states, c_nms], na.rm=T)
	tot_wt_rp_correct = tot_wt_rp * mod$r_to_dsamp[["pool"]][states, "wt"]
	min_wt_rp = min(tot_wt_rp_correct)
	r_state_equalize = min_wt_rp/tot_wt_rp_correct
	all_mats_s = list()
	all_kms = list()
	for(mode in modes) {
		intervs = mod$ehs[mod$ehs$mode == mode,]
		message("processing mode ", mode, " n = ", nrow(intervs))

		for (st in states){
			message("get mat for ", st)
			ratio = r_state_equalize[st]
			spat = get_cmp_mats(mod, state = st, intervs = intervs, dnmt_dko_mod="compare", window=1000, smooth_w = 0)
			
			st_mat = cbind(spat$dnmt*ratio, spat$wt*ratio, spat$tet*ratio)	
			colnames(st_mat) = paste(rep(st,e=wn*3), rep(c("dnmt", "wt", "tet"),e = wn), rep(1:wn, t=3), sep="_")

            if (filter_low_signal){
                pre_f = nrow(st_mat)
                f = rowSums(st_mat) >= 10
                st_mat = st_mat[f,]
                intervs = intervs[f,]
                cli::cli_alert(sprintf("%s loci filtered due to low signal.", pre_f - nrow(st_mat)))
            }

            cli::cli_alert(sprintf("st_mat dim: %s x %s, intervs nrow: %s", nrow(st_mat), ncol(st_mat), nrow(intervs)))
			
			message("will run km")
			k_km = max(min(max_k, floor(nrow(st_mat)/min_n_per_km)),4)
			km = tglkmeans::TGL_kmeans(log2(0.25+st_mat), k=k_km, id_column=F, seed = km_seed)
			rownames(st_mat) = apply(intervs[,1:3], 1, function(x) {x=trimws(x); sprintf("%s:%s_%s", x[1],x[2],x[3])})
			names(km$cluster) = rownames(st_mat)
			mat_s = reord_smooth_on_km(st_mat, km=km, w_size=20, should_hc = T)

			all_mats_s[[st]][[mode]] = mat_s
			all_kms[[st]][[mode]] = km
		}
	}
	return(list(all_mats = all_mats_s, all_kms = all_kms))
}

order_rows_by_weighted_average <- function(mat_subset){
    
    max_per_row = apply(mat_subset, 1, max)
    strong_cl_mat = mat_subset
    
    weighted_avg = as.vector(t(apply(strong_cl_mat, 1, function(row) {
        normalized_row <- (row - min(row)) / (max(row) - min(row))
        wa = sum(normalized_row * seq(101)) / sum(normalized_row)
        return(wa)
    })))
    
    names(weighted_avg) = rownames(strong_cl_mat)    
    strong_cl_order = rownames(strong_cl_mat)[order(weighted_avg)]
    rownames_order = strong_cl_order

    return(rownames_order)
}

order_rows_by_max <- function(mat_subset){
    
    max_per_row = sort(apply(mat_subset, 1, max))
    rownames_order = names(max_per_row)

    return(rownames_order)
}

plot_spat_kms = function(this_mode, st, mats, kms, mode_md, prefix, 
						 order_by = "max", win_n = 101, signal_eps = 0.25)
{
    
    mat = spat_mats[[st]][[this_mode]]
    km = spat_kms[[st]][[this_mode]]

    centers_mat = km$centers
    rownames(centers_mat) = paste0("raw_C", seq(nrow(centers_mat)))

    prefixes = gsub("_[0-9]+", "", colnames(centers_mat))
    max_per_prefix = t(tgs_matrix_tapply(centers_mat, prefixes, max))
    chosen_prefix = names(which.max(colSums(max_per_prefix)))

    prefix_mat = centers_mat[,grepl(chosen_prefix, colnames(centers_mat))]

    if (order_by == "max"){
        raw_C_order <- order_rows_by_max(prefix_mat) 
    } else if (order_by == "weighted_average") {
        raw_C_order <- order_rows_by_weighted_average(prefix_mat)
    } else {cli::cli_abort(message = "order_by must be either 'max' or 'weighted_average'.")}
    
    raw_C_order = as.numeric(gsub("raw_C", "", raw_C_order))
	C_order = paste0("C", rev(seq(length(raw_C_order))))
    cluster_dict = setNames(C_order, raw_C_order)
    mode_md[[this_mode]]$clusters_renaming <<- cluster_dict
    
	temp = km$cluster[km$cluster %in% raw_C_order]
    
    cli::cli_alert(sprintf("%s: %s loci.", this_mode, length(temp)))

    km_cluster = cluster_dict[as.character(temp)] # Translate to new cluster names
	names(km_cluster) = names(temp)
    cluster_factor = factor(km_cluster, levels = C_order)
    cluster_factor = sort(cluster_factor)

    print(sort(table(cluster_factor)))
    
    mat_reordered = mat[names(cluster_factor),]
    shades = colorRampPalette(c("white", "deepskyblue", "deepskyblue4", "black"))(1000)

    
    # Plot spatial matrix
    # 2 versions: one with cluster numbers, one without
    for (plot_version in c("num", "nonnum")){

        show_fig(sprintf("%s%s.%s_%s.png", prefix, this_mode, st, plot_version), height = 1200, width = 600, draw = {
        image(log2(signal_eps+t(mat_reordered)), col=shades, xaxt='n', yaxt='n', zlim=c(-2, 4))
        if (plot_version == "num"){
            mtext(text = gsub("C", "", as.character(C_order)), 
                at=tapply(seq(0,1, l=nrow(mat_reordered)), cluster_factor, mean), 
                side=2, las=2)
        }

        for(i in C_order) {
            abline(h=max(seq(0,1, l=nrow(mat_reordered))[cluster_factor==i]))
        }
        for(i in 1:floor(ncol(mat_reordered)/win_n)) {
            abline(v = i*win_n/ncol(mat_reordered))
        }
        })

    }
        
    return(mode_md)

}

smooth_trend <- function(df_subset, nm, win, smooth_w){
    
    trend = tapply(df_subset[[nm]], df_subset$d_center, function(x) {mean(x, na.rm = TRUE)})
   
    # Fill in missing data points
    missing_dp = setdiff(as.character(seq(-win, win)), names(trend))
    trend[missing_dp] = NA
    trend = trend[order(as.numeric(names(trend)))]
    i_start = which(names(trend) == -win)
    i_end = which(names(trend) == win)

    trend_sm = frollmean(trend, n = smooth_w, na.rm = T)
    trend_sm = trend_sm[i_start:i_end]
    names(trend_sm) = seq(-win, win)
        
    return(trend_sm)
}

extract_meth <- function(intervs, tr_exps, win, smooth_w = NULL){
    
    region_extend = ceiling((smooth_w*20)/2)
    center = intervs
    center$start = center$center
    center$end = center$center + 1
    
        
    gvtrack.create("d_center", center, "distance")
    region = center
    region$start = region$start - (win + region_extend)
    region$end = region$end + (win + region_extend)
            
    spat = gextract(c("d_center"),
                    intervals = region,
                    iterator=20, colnames=c("d_center"))

    spat$d_center_bin = floor(spat$d_center/20) 
    meth = gextract_meth(tracks = tr_exps, names = names(tr_exps), intervals = region, iterator = 20)
    
    spat_meth = gintervals.neighbors1(spat, meth, maxdist=0, na.if.notfound = T)
    spat_meth = spat_meth[,c("chrom", "start", "end", "intervalID", "d_center", "d_center_bin", names(tr_exps))]
    spat_meth = spat_meth %>% filter(abs(d_center) <= win)
    
    return(spat_meth)
}

plot_meth_single_origin <- function(mod, mode, km, st, meth_tracks, win = 1e+3, smooth_w = NULL, panel = NULL){

    mode_intervs = mod$ehs[mod$ehs$mode == mode,]
    mode_f = data.frame(items = rownames(spat_mats[[st]][[mode]])) %>%
        separate(items, into = c("chrom", "start_end"), sep = ":") %>%
        separate(start_end, into = c("start", "end"), sep = "_") %>%
        mutate(start = as.numeric(start), end = as.numeric(end))

    mode_intervs = mode_intervs %>% right_join(mode_f)
    message(mode)

    mode_intervs = mode_intervs[,c("chrom", "start", "end", "strand", "center")]
    mode_intervs = mode_intervs %>% mutate(intervalID = row_number(), 
                                           cluster = as.character(km[[st]][[mode]]$cluster))

    spat = extract_meth(intervs = mode_intervs, tr_exps = meth_tracks, win = win, smooth_w = smooth_w)

    spat = spat %>% left_join(mode_intervs %>% select(intervalID, cluster)) 
    n_cluster = mode_intervs %>% group_by(cluster) %>% summarize(n = n())

    trends_list = list()
    for (nm in names(meth_tracks)){
        trend_l <- lapply(split(spat, spat$cluster), function(x) {smooth_trend(x, nm = nm, win, smooth_w)})
        trend_df <- as.data.frame(do.call(cbind, trend_l))                                        
        trend_df_long = trend_df %>% rownames_to_column("d_center") %>% mutate(d_center = as.numeric(d_center)) %>% pivot_longer(-d_center, names_to = "cluster", values_to = "meth") %>%
            left_join(n_cluster) %>% mutate(facet = sprintf("C%s: n = %s", cluster, n)) %>% arrange(as.numeric(cluster)) %>% mutate(origin = nm)
    trends_list[[nm]] = trend_df_long
    
    }

    plot_df = do.call(rbind, trends_list)

    cl_plots = list()
    for (i in n_cluster$cluster){

        p <- plot_df %>% filter(cluster == i) %>% 
            ggplot(aes(x = d_center, y = meth, color = origin, linetype = origin)) + geom_line() + ylim(-0.1, 1) +
                theme_minimal() + theme(legend.position = "None", plot.title = element_text(size = 14, face = "bold")) + 
                xlab("") + ylab("") + 
                geom_hline(yintercept = 0, color = "blue", lty = 2) + 
                scale_color_manual(values = c("wt" = "black", "tet" = "red", "tet_dko" = "darkorange", "tet1cd" = "chocolate3")) +
                scale_linetype_manual(values = c("wt" = "dashed", "tet" = "solid", "tet_dko" = "solid", "tet1cd" = "solid"))

        cl_plots[[i]] = p

    }
    
    return(cl_plots)
}

plot_acc_single_origin <- function(df, origin_to_highlight, custom_ymax){
    ymin = -2
    ymax = custom_ymax
    df_f = df %>% filter(origin %in% c("wt", origin_to_highlight))
    p <- ggplot(df_f, aes(d_center, acc, color = origin, linetype = origin)) + geom_line() +
        scale_color_manual(values = c("wt" = "black", "dnmt" = "blue", "tet" = "red")) +
        scale_linetype_manual(values = c("wt" = "dashed", "dnmt" = "solid", "tet" = "solid")) + 
        theme_minimal() + theme(legend.position = "None") + xlab("") + ylab("") + ylim(ymin, ymax)
    return(p)
}

plot_acc_trend <- function(st, this_mode, fig_format = "pdf"){

    win_n = 101
    offs = seq(-20*(win_n-1)/2, 20*(win_n-1)/2,20)

    km_centers = spat_kms[[st]][[this_mode]]$centers
    rownames(km_centers) = 1:max(spat_kms[[st]][[this_mode]]$cluster)

    km_centers_dnmt = km_centers[,1:win_n]
    km_centers_wt = km_centers[,(win_n+1):(2*win_n)]
    km_centers_tet = km_centers[,(2*win_n+1):(3*win_n)]

    cluster_list = list()
    for (i in 1:max(spat_kms[[st]][[this_mode]]$cluster)){
        cluster_long = data.frame("d_center" = seq(-20*(win_n-1)/2, 20*(win_n-1)/2,20),
                             "wt" = km_centers_wt[i,],
                             "dnmt" = km_centers_dnmt[i,],
                             "tet" = km_centers_tet[i,]) %>% 
            pivot_longer(-d_center, names_to = "origin", values_to = "acc") %>% 
            mutate(cluster_ind = i)

        cluster_list[[i]] = cluster_long
    }

	plot_df = do.call(rbind, cluster_list)
    plot_df = plot_df %>%
    	mutate(cluster_ind = factor(cluster_ind, levels = sort(unique(plot_df$cluster_ind))))    
	if (this_mode %in% c("tss_hcg", "ins_sctcf", "icg_senh")){
        custom_ymax = 5} else if (this_mode %in% c("hcg_senh", "icg_wenh", "tss_icg", "tss_lcg")){
        custom_ymax = 2} else if (this_mode == "tss_icg"){
		custom_ymax = 3} else {cli::cli_abort(message = "Unrecognized mode name.")}

    cluster_dict = mode_md[[this_mode]]$clusters_renaming
    for (i in names(cluster_dict)){
        j = cluster_dict[[i]]
        plot_df_subset = plot_df %>% filter(cluster_ind == i)

        p_dnmt = plot_acc_single_origin(plot_df_subset, "dnmt", custom_ymax)
        p_tet = plot_acc_single_origin(plot_df_subset, "tet", custom_ymax)

        p_both = ggarrange(p_dnmt, p_tet)
        p_both <- annotate_figure(p_both, top = text_grob(paste("Cluster", gsub("C", "", j)), 
                        face = "bold", size = 14))
    
        if (!dir.exists(sprintf("./figs/Spatials/acc_trends/%s/", this_mode))){
            dir.create(out_path(sprintf("./figs/Spatials/acc_trends/%s/", this_mode)))
        }
        
        # If fig format is "png", save as png, otherwise save as pdf
        if (fig_format == "png"){
            show_fig(sprintf("./figs/Spatials/acc_trends/%s/%s.png", this_mode, j), p_both, show = FALSE, width = 5, height = 2)
        } else {
            show_fig(sprintf("./figs/Spatials/acc_trends/%s/%s.pdf", this_mode, j), p_both, show = FALSE, width = 5, height = 2)
        }

    }
}

plot_meth_trend <- function(st, this_mode, fig_format = "pdf"){
    
    cluster_dict = mode_md[[this_mode]]$clusters_renaming
    cl_plots_1 = plot_meth_single_origin(mod, this_mode, km = spat_kms, st, 
                    meth_track = c("wt" = "CHIM.comb_wt_E8_5",
                                   "tet" = "CHIM.comb_tet_tko_E8_5"),
                    smooth_w = 20, panel = 1)

    cl_plots_2 = plot_meth_single_origin(mod, this_mode, km = spat_kms, st, 
                    meth_track = c("wt" = "embdnmt.HTNC",
                                   "tet_dko" = "embdnmt.DKO23S1",
                                   "tet1cd" = "embdnmt.Tet1CDS1"), 
                    smooth_w = 20, panel = 2)
    
    for (i in names(cluster_dict)){
        
        j = cluster_dict[[i]]
        p1 = cl_plots_1[[i]]
        p2 = cl_plots_2[[i]]
        p_both = ggarrange(p1, p2)
        p_both <- annotate_figure(p_both, top = text_grob(paste("Cluster", gsub("C", "", j)), 
                                  face = "bold", size = 14))
        
        # If fig format is "png", save as png, otherwise save as pdf
        if (fig_format == "png"){
            show_fig(sprintf("./figs/Spatials/meth_trends/%s/%s.png", this_mode, j), p_both, show = FALSE, width = 5, height = 2)
        } else {
            show_fig(sprintf("./figs/Spatials/meth_trends/%s/%s.pdf", this_mode, j), p_both, show = FALSE, width = 5, height = 2)
        }    

    }
}

# ---- used by Fig4.ipynb ----

plot_spat_kms_named <- function (this_mode, st, mats, kms, mode_md, prefix, order_by = "max", 
    win_n = 101, signal_eps = 0.25) 
{
    mat = spat_mats[[st]][[this_mode]]
    km = spat_kms[[st]][[this_mode]]
    centers_mat = km$centers
    rownames(centers_mat) = paste0("raw_C", seq(nrow(centers_mat)))
    prefixes = gsub("_[0-9]+", "", colnames(centers_mat))
    max_per_prefix = t(tgs_matrix_tapply(centers_mat, prefixes, 
        max))
    chosen_prefix = names(which.max(colSums(max_per_prefix)))
    prefix_mat = centers_mat[, grepl(chosen_prefix, colnames(centers_mat))]
    if (order_by == "max") {
        raw_C_order <- order_rows_by_max(prefix_mat)
    }
    else if (order_by == "weighted_average") {
        raw_C_order <- order_rows_by_weighted_average(prefix_mat)
    }
    else {
        cli::cli_abort(message = "order_by must be either 'max' or 'weighted_average'.")
    }
    raw_C_order = as.numeric(gsub("raw_C", "", raw_C_order))
    C_order = paste0("C", rev(seq(length(raw_C_order))))
    cluster_dict = setNames(C_order, raw_C_order)
    ###
    mode_md[[this_mode]]$clusters_renaming <<- cluster_dict
    ###
    temp = km$cluster[km$cluster %in% raw_C_order]
    cli::cli_alert(sprintf("%s: %s loci.", this_mode, length(temp)))
    km_cluster = cluster_dict[as.character(temp)]
    names(km_cluster) = names(temp)
    cluster_factor = factor(km_cluster, levels = C_order)
    cluster_factor = sort(cluster_factor)
    mat_reordered = mat[names(cluster_factor), ]
    ###
    mode_md[[this_mode]]$intervs_order <<- rownames(mat_reordered)
    mode_md[[this_mode]]$cluster_factor <<- cluster_factor
    ###
    shades = colorRampPalette(c("white", "deepskyblue", "deepskyblue4", 
        "black"))(1000)
    for (plot_version in c("num", "nonnum")) {
        show_fig(sprintf("%s%s.%s_%s.png", prefix, this_mode, st, 
            plot_version), height = 1200, width = 600, draw = {
        image(log2(signal_eps + t(mat_reordered)), col = shades, 
            xaxt = "n", yaxt = "n", zlim = c(-2, 4))
        if (plot_version == "num") {
            mtext(text = gsub("C", "", as.character(C_order)), 
                at = tapply(seq(0, 1, l = nrow(mat_reordered)), 
                  cluster_factor, mean), side = 2, las = 2)
        }
        for (i in C_order) {
            abline(h = max(seq(0, 1, l = nrow(mat_reordered))[cluster_factor == i]))
        }
        for (i in 1:floor(ncol(mat_reordered)/win_n)) {
            abline(v = i * win_n/ncol(mat_reordered))
        }
        })
    }
    return(mode_md)
}

# Factor that scales each track to the smallest library of the set
scale_to_min <- function(tracks){
    totals <- setNames(vapply(tracks, \(tr) unname(gsummary(tr)["Sum"]), 1), tracks)
    min(totals, na.rm = TRUE) / totals
}

# Sum each track in 280bp windows at 20bp steps and write a UCSC wig (scale, if given, multiplies track tr by scale[[tr]])
write_gb_wigs <- function(tracks_to_extract, scale = NULL, span = 20){
    names(track_colors) = tracks_to_extract
    for (tr in tracks_to_extract){
        gvtrack.create(paste0("v_",tr), tr, "sum")    
        gvtrack.iterator(paste0("v_",tr), sshift = -140, eshift = 140)
        expr = if (is.null(scale)) paste0("v_",tr) else paste0("v_",tr, " * ", scale[[tr]])
        data <- gextract(expr, gintervals.all() %>% filter(chrom != "chrM"), iterator = span, colnames = "score")
        dataf = data %>% left_join(gintervals.all() %>% dplyr::rename(chrstart=start, chrend=end), by="chrom") %>% filter(end+span <= chrend)
        fwrite_ucsc(dataf, file = file.path(ucsc_dir, paste0(tr, ".wig")), name = tr, description = tr, type="wig", span = span, color = track_colors[[tr]])
    }
    for (tr in tracks_to_extract){
        cat(sprintf("wigToBigWig %s/%s.wig %s %s/%s.bw\n", ucsc_dir, tr, mm10_chrom_sizes, ucsc_dir, tr))
    }
}
