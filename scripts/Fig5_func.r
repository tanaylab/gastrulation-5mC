
library("anndata")
library("tglkmeans")
library(tibble)
library(tidyr)

load_expression = function(mod, recomp=F)
{
	if(file.exists("data/output/rna_stat.RDS") & !recomp) {
		mod$rna_stat = readRDS(in_path("data/output/rna_stat.RDS"))
		return(mod)
	}

    message("will read cells")
    cells = read_h5ad("./data/MC_models/RNA_all.h5ad")	
    message("done read cells")
    umis = t(as.matrix(cells$X))
    message("done cells 2 mat")
    n_umi = colSums(umis)

    cells$obs$cell_name = rownames(cells$obs)
    colnames(umis) = cells$obs$cell_name

    cell_host = cells$obs$KO_host
    names(cell_host) = cells$obs$cell_name
    cell_geno = as.character(cells$obs$genotype)
    names(cell_geno) = cells$obs$cell_name
    cell_genoko = paste(cell_geno, cell_host, sep="_")
    names(cell_genoko) = cells$obs$cell_name
    cell_embko = paste(cell_host, cells$obs$embryo, sep="_")
    names(cell_embko) = cells$obs$cell_name
    cell_emb = as.character(cells$obs$embryo)
    names(cell_emb) = cells$obs$cell_name

    cell_annots = read.csv("data/ref/cell_annotation.csv", stringsAsFactors=F)
    type_layer = read.table("data/ref/annot_key.txt", stringsAsFactors=F,h=T, sep="\t")

	cell_mc_type = cell_annots$cell_type
	names(cell_mc_type) = cell_annots$cell

	type_to_layer = type_layer$layer
	names(type_to_layer) = type_layer$type

	cell_layer = type_to_layer[cell_mc_type]
	names(cell_layer) = cell_annots$cell

	layer_umis = list()
	layers = unique(cell_layer)
	layers = setdiff(layers, "NA")
	layers = layers[!is.na(layers)]
	genotypes = unique(cell_geno[names(cell_layer)])
	genotypes = genotypes[!is.na(genotypes)]

	margs_ko = list()
	margs_host = list()
	for(geno in genotypes) {
		marg_ko= NULL
		marg_host= NULL
		act_layer = NULL
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
			if(!is.null(ds_umis)) {
				layer_umis[[sprintf("%s.%s", geno, layer)]] = ds_umis
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

	mko_ds = downsample_matrix(mko, 250000)
	rownames(mko_ds) = rownames(mko)
	colnames(mko_ds) = colnames(mko)
	mhost_ds = downsample_matrix(mhost, 250000)
	rownames(mhost_ds) = rownames(mhost)
	colnames(mhost_ds) = colnames(mhost)

	mko_e = t(t(mko_ds)/colSums(mko_ds))
	mhost_e = t(t(mhost_ds)/colSums(mhost_ds))

	mko_l = log2(2e-5+mko_e)
	mhost_l = log2(2e-5+mhost_e)
	mko_lf = mko_l - apply(mko_l,1,median)
	mhost_lf = mhost_l - apply(mhost_l,1,median)

	mdlt = mko_l - mhost_l

	mko_fe = t(t(mko)/colSums(mko))
	mhost_fe = t(t(mhost)/colSums(mhost))
	mdlt_full = log2(1e-5++mko_fe)-log2(1e-5+mhost_fe)

	gsetroot(mm10_trackdb)
	tss = gintervals.load("intervs.global.tss")
	gene_chrom = as.character(tss$chrom[ ! duplicated(tss$geneSymbol)])
	names(gene_chrom) = tss$geneSymbol[! duplicated(tss$geneSymbol)]

	g_cg = gextract(c("seq.CG_500_mean_new"), intervals=tss, iterator=tss, colnames=c("cg"))
	gene_cg = g_cg$cg
	names(gene_cg) = tss$geneSymbol

	tab_weber = read.csv("data/weber_genes.csv",stringsAsFactors=F)
	weber_germ = intersect(rownames(umis), tab_weber$gene[tab_weber$class=="upregulated_germline"])
	weber_other = intersect(rownames(umis), tab_weber$gene[tab_weber$class!="upregulated_germline"])
	imprints = c("Igf2", "Dlk1", "Meg3", grep("Peg", rownames(umis),v=T), "Dnmt3a", "Dnmt3b", "Dnmt1", "Mest", "H19", "Impact","Airn","Kcnq1ot1", "Snrpn","Nnat","Rian","Inpp5f","Plagl1","Nap1l5","Gnas","Cdkn1c","Phlda2")

	g_not_6x = intersect(names(gene_chrom), rownames(mko))
	g_not_6x = setdiff(g_not_6x, names(which(gene_chrom == "chr6")))
	g_not_6x = setdiff(g_not_6x, names(which(gene_chrom == "chrX")))
	g_tko = g_not_6x
	g_not_1x = intersect(names(gene_chrom), rownames(mko))
	g_not_1x = setdiff(g_not_1x, names(which(gene_chrom == "chr1")))
	g_not_1x = setdiff(g_not_1x, names(which(gene_chrom == "chrX")))
	g_dko = g_not_1x

	gboth = intersect(g_tko, g_dko)

	ty_stats = list()
	for(ty in c("Ecto","EmbMeso", "ExMeso", "Endothel", "Ery", "Endoderm")) {
		tko_nm = sprintf("DNMT_TKO.%s", ty)
		dko_nm = sprintf("DNMT_DKO.%s", ty)

		cmptko = cbind(mko[gboth,tko_nm],mhost[gboth,tko_nm])
		tot_tko = colSums(cmptko)
		cmpdko = cbind(mko[gboth,dko_nm],mhost[gboth,dko_nm])
		tot_dko = colSums(cmpdko)
		pv3 = apply(cmptko, 1, function(x) return(chisq.test(rbind(x, tot_tko-x))$p.value))
		pv2 = apply(cmpdko, 1, function(x) return(chisq.test(rbind(x, tot_dko-x))$p.value))
		pv3[is.nan(pv3)] = 1
		pv2[is.nan(pv2)] = 1
		e3 = log2(2e-6+t(t(cmptko)/tot_tko))
		e2 = log2(2e-6+t(t(cmpdko)/tot_dko))
		fc3 = e3[,1]-e3[,2]
		fc2 = e2[,1]-e3[,2]
		stats = cbind(cbind(fc3, fc2), cbind(pv3, pv2))
		colnames(stats) = paste(rep(ty,4), c("lfTKO", "lfDKO", "pvTKO", "pvDKO"), sep=".")

		sig = pv3<1e-6 & pv2 < 1e-6
		message("at ", ty, " sig ", sum(sig))

		ty_stats[[ty]] = stats
	}

	de = do.call("cbind", ty_stats)
	de = as.data.frame(de)

	rs = list(gboth=gboth, ty_stats = ty_stats, weber_germ=weber_germ, weber_other=weber_other, imprints = imprints, gene_cg = gene_cg, mko = mko, mhost=mhost, mdlt=mdlt, mdlt_full=mdlt_full)
	mod$rna_stat = rs
	saveRDS(rs, out_path("data/output/rna_stat.RDS"))
	return(mod)
}

calc_chisq_pval <- function(umi_df_long, UMI_thresh = 10, num_proc = 40){

    doMC::registerDoMC(num_proc)
    chisq_results <- umi_df_long %>% 
		mutate(under_rep = wt + dnmt < UMI_thresh) %>%
        plyr::adply(1, function(x) {
            chisq <- chisq.test(
                matrix(c(x$dnmt, x$tot_min_dnmt, x$wt, x$tot_min_wt), byrow = TRUE, nrow = 2)
            )
            tibble(pval = chisq$p.value, statistic = chisq$statistic)
        }, .parallel = TRUE)

    return(chisq_results)
}

chisq_raw_atac <- function(raw_atac){
    raw_atac_f = raw_atac[,1:2]
    colnames(raw_atac_f) = c("wt", "dnmt")
    raw_long = raw_atac_f %>% as.data.frame() %>% mutate(tot_min_wt = sum(wt) - wt, tot_min_dnmt = sum(dnmt) - dnmt)
    chisq = calc_chisq_pval(raw_long)
    
    chisq = chisq %>% mutate(pval = ifelse(under_rep | is.nan(pval), NA, pval),
                             statistic = ifelse(under_rep | is.nan(pval), NA, pval))
    
    return(chisq)
}

plt_exp_atac_tss = function(mod){
	f_tss = abs(mod$npeaks$d_tss) < 1e+3

    rs = mod$rna_stat
    mdlt = as.data.frame(rs$mdlt_full)
    mdlt$geneSymbol = rownames(mdlt)

    tss = gintervals.load("intervs.global.tss")
    tss_cov = tss[tss$geneSymbol %in% rownames(rs$mdlt),c(1:4,9)]
    tss_cov = tss_cov %>% left_join(mdlt)

    lp = mod$lin_peaks

    u_rnp = lp$rnp_ds[f_tss,]
    e_rnp  = t(t(u_rnp)/colSums(u_rnp))
    le_rnp = log2(1e-6+e_rnp)

    u_acp = lp$acp_ds[f_tss,]
    e_acp  = t(t(u_acp)/colSums(u_acp))
    le_acp = log2(1e-6+e_acp)

    u_meso = lp$meso_ds[f_tss,]
    e_meso  = t(t(u_meso)/colSums(u_meso))
    le_meso = log2(1e-6+e_meso)

    tss_peaks = mod$npeaks[f_tss,]

    # generate a matrix on TSS with dnmt-host rna, dnmt-wt atac
    tss_merge = gintervals.neighbors(tss_peaks, tss_cov)

    f_hit = tss_merge$dist==0
    tss_merge = tss_merge[f_hit,]
    le_rnp = le_rnp[f_hit,]
    le_acp = le_acp[f_hit,]
    le_meso = le_meso[f_hit,]
    tss_peaks = tss_peaks[f_hit,]
    f_icg = tss_peaks$CG<0.04 & tss_peaks$CG>0.02
    f_lcg = tss_peaks$CG<0.02
    f_hcg = tss_peaks$CG>0.04

    # ATAC diff
    a_stat = data.frame(a_rnp = le_rnp[,2]-le_rnp[,1], 
                        a_acp = le_acp[,2]-le_acp[,1], 
                        a_meso = le_meso[,2]-le_meso[,1])
    
	a_mean_diff = apply(a_stat, 1, mean)					
    # Mean ATAC
    abs_a = (le_rnp[,1]+le_acp[,1]+le_meso[,1])/3

    # RNA diff
    r_stat = data.frame(
            e_rnp = tss_merge$DNMT_TKO.Ecto+tss_merge$DNMT_DKO.Ecto,
            e_acp = tss_merge$DNMT_TKO.ExMeso+tss_merge$DNMT_DKO.ExMeso,
            e_meso = tss_merge$DNMT_TKO.EmbMeso+tss_merge$DNMT_DKO.EmbMeso)
    gnms = tss_merge$geneSymbol

    # Remove duplicate genes, choose the one with the highest a_stat
    is_to_remove = c()
    dup_genes = gnms[duplicated(gnms)]
    for (dg in dup_genes){
        is = which(gnms == dg)
        max_i = is[which.max(a_mean_diff[is])]
        
        is_to_remove = c(is_to_remove, setdiff(is, max_i))
    }

    a_stat = a_stat[-is_to_remove,]
    r_stat = r_stat[-is_to_remove,]
    abs_a = abs_a[-is_to_remove]
    gnms = gnms[-is_to_remove]
	tss_merge = tss_merge[-is_to_remove,]

    # Remove imprinted genes to be plotted separately
    is_imprinted = which(gnms %in% setdiff(mod$rna_stat$imprints, c("Dnmt1", "Dnmt3a", "Dnmt3b")))
    a_stat_imprinted = a_stat[is_imprinted,]
    r_stat_imprinted = r_stat[is_imprinted,]
    abs_a_imprinted = abs_a[is_imprinted]
    gnms_imprinted = gnms[is_imprinted]

    is_imprinted = which(gnms %in% mod$rna_stat$imprints)
    a_stat = a_stat[-is_imprinted,]
    r_stat = r_stat[-is_imprinted,]
    abs_a = abs_a[-is_imprinted]
    gnms = gnms[-is_imprinted]
	tss_merge = tss_merge[-is_imprinted,]
    
	#########################
	# STATISTICAL TEST - ATAC
	# Filter ATAC p-values using the same filtering regime as a_stat

	if (!file.exists("./data/output/atac_chisq.RDS")){
		rnp_chisq = chisq_raw_atac(lp$rnp_ds)
		meso_chisq = chisq_raw_atac(lp$meso_ds)
		acp_chisq = chisq_raw_atac(lp$acp_ds)

		atac_chisq = list("rnp_chisq" = rnp_chisq, "meso_chisq" = meso_chisq, "acp_chisq" = acp_chisq)
		saveRDS(atac_chisq, out_path("./data/output/atac_chisq.RDS"))
	} else {
		atac_chisq = readRDS(in_path("./data/output/atac_chisq.RDS"))
	}

	filter_atac_chisq = function(chisq, f_tss, f_hit, is_to_remove, is_imprinted){
		chisq = chisq[f_tss,]
		chisq = chisq[f_hit,]
		chisq = chisq[-is_to_remove,]
		chisq = chisq[-is_imprinted,]

		return(chisq)
	}

	rnp_chisq = filter_atac_chisq(atac_chisq$rnp_chisq, f_tss, f_hit, is_to_remove, is_imprinted)
	acp_chisq = filter_atac_chisq(atac_chisq$acp_chisq, f_tss, f_hit, is_to_remove, is_imprinted)
	meso_chisq = filter_atac_chisq(atac_chisq$meso_chisq, f_tss, f_hit, is_to_remove, is_imprinted)

	atac_pvals = data.frame(rnp = rnp_chisq$pval, acp = acp_chisq$pval, meso = meso_chisq$pval)
	
	# Run FDR correction on the entire dataframe
	pvals_vector <- as.vector(as.matrix(atac_pvals))
	pval_fdr_vector = rep(NA, length(pvals_vector)) # Initialize with NAs
	pval_fdr_vector[!is.na(pvals_vector)] <- p.adjust(pvals_vector[!is.na(pvals_vector)], method = "fdr")

	atac_qvals <- matrix(pval_fdr_vector, nrow = nrow(atac_pvals), ncol = ncol(atac_pvals))                                        
	rownames(atac_qvals) <- rownames(atac_pvals)

	# ATAC deltas with p-values > 0.01 are set to 0
	a_stat[atac_qvals > 0.01] = 0
    
	#########################
	#########################
	# STATISTICAL TEST - RNA
	RNA_chisq = readRDS(in_path("./data/output/chisq_screen.RDS"))
	rna_pvals = RNA_chisq %>% filter(layer %in% c("EmbMeso", "ExMeso", "Ecto"), geno == "DNMT_TKO")
	rna_pvals$qval = p.adjust(rna_pvals$pval, method = "fdr")
	rna_qvals = rna_pvals %>% select(gene, layer, qval) %>% pivot_wider(id_cols = gene, names_from = layer, values_from = qval, values_fill = 1) %>% column_to_rownames("gene")

	rna_qvals = rna_qvals[gnms,c("Ecto", "ExMeso", "EmbMeso")] %>% mutate_all(~ replace_na(., 1)) %>% as.matrix()

	# RNA deltas with p-values > 0.01 are set to 0
	r_stat[rna_qvals > 0.01] = 0
	#########################

    # Filter for strong diff in ATAC/RNA
    fa = apply(a_stat, 1, max) > 0.8
    fr = apply(r_stat, 1, max) > 1
	
	# message(sprintf("fa: %s genes.",sum(fa)))
	# message(sprintf("fr: %s genes.",sum(fr)))

	atac_qvals[a_stat <= 0.8] = NA

    fr_d = apply(r_stat, 1, min) < -1

	bshades = colorRampPalette(c("white", "lightblue1", "lightblue2", "deepskyblue2", "deepskyblue4", "black"))(1000)
	tshades = colorRampPalette(c("white", "cadetblue1", "cadetblue2", "cadetblue3", "cadetblue4"))(1000)
	gshades = colorRampPalette(c("white", "grey95", "darkolivegreen3", "darkolivegreen4", "darkolivegreen", "black"))(1000)

    for(lab in c("11","01","10")) {
        if (lab == "11") { 
            f = fa & fr # Strong diff both  
        } else if(lab == "01") { 
            f = !fa & fr # Strong diff RNA  
        } else if(lab == "10") {
            f = fa & !fr # Strong diff ATAC
        } else {
            f = fr_d # Down
        }

        km = tglkmeans::TGL_kmeans(cbind(a_stat, r_stat)[f,], k=10, id_column=F, seed=19)
		kord = order(km$cluster + abs_a[f]/1000)
		# Save kord as RDS
		saveRDS(kord, out_path(sprintf("./data/output/5C_kord_%s.RDS", lab)))

        gaps_row = which(!duplicated(sort(km$cluster)))
        gaps_row = gaps_row[2:length(gaps_row)]
        gaps_row = gaps_row - 1
        
		# Write out the cluster assignments
		tss_coord = tss_merge[,1:3]
		data_out = cbind(tss_coord[f,][kord,], gnms[f][kord], km$cluster[kord])
		colnames(data_out) = c("chrom", "start", "end", "gene", "cluster")
        write.table(data_out, out_path(sprintf("data/summary/tss_km_%s.txt",lab)), quote=F, sep="\t")
        
        message(sprintf("%s: %s genes.", lab, length(gnms[f][kord])))

		gns_to_highlight = read.csv("./data/ref/5C_highlight.csv") %>% pull(1)
		labels_row <- gnms[f][kord]
		labels_row[!grepl(paste(gns_to_highlight, collapse = "|"), labels_row, ignore.case = TRUE)] = ""
        
        show_fig(sprintf("figs/Fig5/final_heatmaps/tss_rna_%s.no_gaps.pdf", lab), draw = pheatmap::pheatmap(r_stat[f,][kord,], 
                        cluster_rows=F, cluster_cols=F, 
                        w=4, h=sum(f)/100, fontsize = 15, 
                        col=tshades, labels_row = labels_row,
                        legend = F,
                        show_colnames=F, breaks = seq(0,3,l=1001)))
        
        show_fig(sprintf("figs/Fig5/final_heatmaps/tss_atac_%s.no_gaps.pdf", lab), draw = pheatmap::pheatmap(a_stat[f,][kord,], 
                        cluster_rows=F, cluster_cols=F, 
                        w=4, h=sum(f)/100, fontsize = 15, 
                        col=bshades, labels_row = labels_row,
                        legend = F,
                        show_colnames=F, breaks = seq(0,3,l=1001)))
        
        show_fig(sprintf("figs/Fig5/final_heatmaps/tss_abs_a_%s.no_gaps.pdf", lab), draw = pheatmap::pheatmap(as.matrix(abs_a[f][kord]), 
                        cluster_rows=F, cluster_cols=F, 
                        w=1, h=sum(f)/100, fontsize = 15,
                        col=gshades, show_rownames=F,
                        legend = F,
                        show_colnames=F, breaks = seq(-20,-13,l=1001)))
            
    return(mod)
    }

}

# ---- used by Fig5.ipynb ----

generate_marg_df <- function(intervs){

    ### DNMT TKO
    gvtrack.create("TKO", "mut_multi.DNMT_TKO_KO.raw", "sum")
    TKO = gextract("TKO", intervs, iterator = intervs)
    TKO$TKO[is.na(TKO$TKO)] = 0
    ### DNMT DKO
    gvtrack.create("DKO", "mut_multi.DNMT_DKO_KO.raw", "sum")
    DKO = gextract("DKO", intervs, iterator = intervs)
    DKO$DKO[is.na(DKO$DKO)] = 0
    ### TET TKO
    gvtrack.create("TET", "mut_multi.TET_TKO_KO.raw", "sum")
    TET = gextract("TET", intervs, iterator = intervs)
    TET$TET[is.na(TET$TET)] = 0
    gvtrack.create("WT_lg", "mut_multi.WT_lg.raw", "sum")
    gvtrack.create("WT_sm", "mut_multi.WT_sm.raw", "sum")
    WT = gextract("WT_lg + WT_sm", intervs, iterator = intervs, colnames = "WT")
    WT$WT[is.na(WT$WT)] = 0
    
    df_all = cbind(intervs, data.frame("tko" = TKO$TKO,
                                        "dko" = DKO$DKO,
                                        "tet" = TET$TET, 
                                        "wt" = WT$WT))
    
    return(df_all)
}

process_atac <- function(df_all, tr_karyo, r_to_dsamp){

    LTR_proc = df_all
    for (nm in c("tko", "dko", "tet")){
        karyo_r = tr_karyo[[nm]]
        bad_chr = intersect(names(which(karyo_r>1.06)), unique(as.character(LTR_proc$chrom)))
        for(chr in bad_chr) {
            f = LTR_proc$chrom == chr
            r = karyo_r[chr]
            vec = matrix(LTR_proc[f, nm], ncol=1)
            message("correct karyo ", nm, " chr ", chr, " r = ", r, " n = ", round(sum(vec)/r))
            vec_ds = downsample_matrix(vec, target_n=round(sum(vec)/r), seed=28)
            LTR_proc[f,nm] = vec_ds
        }
    }

    for(nm in colnames(r_to_dsamp)) {
        if(!is.null(r_to_dsamp[,nm])) {
            to_ds = matrix(LTR_proc[,nm],ncol=1)
            r_ds = r_to_dsamp[,nm]
            message("into dsamp, ", nm, " targ ", floor(sum(to_ds)*r_ds))
            post_ds = downsample_matrix(to_ds,  target_n = floor(sum(to_ds)*r_ds), seed=19)
            LTR_proc[,nm] = post_ds
        }
    }

    return(LTR_proc)
    
}

plot_largest_families <- function(df_proc){
    p = df_proc %>% filter(name %in%
        (df_proc %>% group_by(name) %>% summarise(name_n = n()) %>% slice_max(name_n, n = 5) %>% pull(name))
    ) %>%
    select(name, wt, tko, dko, tet) %>% pivot_longer(-name, names_to = "origin",) %>%
    mutate(origin = factor(origin, c("tko", "dko", "tet", "wt")), value = log2(0.25+value)) %>%
    ggplot(aes(value, name, fill = origin)) + geom_boxplot(outlier.size = 0.4) + theme_minimal(base_size = 20) + 
        scale_fill_manual(values = c("wt" = "darkolivegreen3", "tko" = "deepskyblue4", "dko" = "deepskyblue2", "tet" = "darkred")) +
        xlab("") + ylab("") +
        theme(legend.position = "none") + 
        xlim(0,10)
    
    return(p)
}

# ---- used by Fig5.ipynb ----

plot_paired_bars <- function(highlight_gn, upper_y){
    df_rect <- df_long %>%
      filter(gene == highlight_gn) %>%
      mutate(
        x_center = if_else(group == "host", -dodge_w/2, dodge_w/2),
        xmin = x_center - bar_w/2,
        xmax = x_center + bar_w/2,
        ymin = lower_y,
        ymax = if_else(value < lower_y, lower_y + min_h, value),
        label_y = pmin(ymax, upper_y) - 0.05
      )

    p <- ggplot(df_rect) +
      geom_hline(yintercept = lower_y, size = 0.3, color = "grey70") +
      geom_rect(
        aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = group),
        color = "grey20",
        size = 0.25
      ) +
      scale_fill_manual(values = c("KO" = "deepskyblue", "host" = "grey50")) +
      scale_x_continuous(limits = c(-0.5, 0.5), expand = c(0, 0)) +
      scale_y_continuous(limits = c(lower_y, upper_y), expand = c(0, 0)) +
      labs(title = paste(highlight_gn, "in TKO"), x = NULL, y = NULL) +
      theme_classic(base_size = 11) +
      theme(
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
          axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
        legend.position = "none"
      )
    
    return(p)
}

plot_tss11 <- function(cols, color){
    plot_df = sc$mdlt_full[,cols]
    colnames(plot_df) = c("TKO", "DKO")
    plot_df = plot_df %>% as.data.frame() %>% rownames_to_column("gene")

    p <- ggplot() + geom_point(data = subset(plot_df, !gene %in% tss11$gene), aes(DKO, TKO), color = "grey30") +
        geom_point(data = subset(plot_df, gene %in% tss11$gene), aes(DKO, TKO), color = color) +
        theme_minimal(base_size = 20) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +
        xlim(-3,3) + ylim(-3,3)
    
    show_fig(sprintf("./figs/Supp/%s_%s.pdf", cols[[1]], cols[[2]]), p, width = 5, height = 5)
    return(invisible(p))
}
