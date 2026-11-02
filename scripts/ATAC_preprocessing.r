library("misha")
library("dplyr")
library("zoo")
library("tglkmeans")
library("misha.ext")
library("anndata")
library("tgstat")

gsetroot(mm10_trackdb)

mark_gs = c("Utf1", "Pou3f1", "Pou5f1", "Dppa2", "Dppa4", "Dppa5a", "Sox2", "Sfrp1", "Hand1", "Igf2", "Lefty1", "Lefty2", "Nodal", "Tdgf1", "Sp5", "Eomes", "Hhex", "Sox1", "Nkx1-2", "Tfap2a", "Meis1", "Mef2c","Alx","Prrx2", "Foxc1","Nanog", "Gsc", "Foxa1", "Foxa2", "Foxj1", "Hhex",
"Mesp1", "Mixl1", "Fgf3", "Aldh1a2", "Crabp1", "Crabp2", "Bmp4", "Msx1", "Kdr", "Tal1", "Runx1", "Foxa2", "Sox17", "Shisa2", "Gata4", "Gata2", "Wnt5a", "Cdkn1c", "Apln", "Cdx1", "Foxc1", "Foxc2", "Twist1", "Tcf15", "Dll1", "Dll3", "Sox2", "Sox3", "Igfbp2", "Sall2", "Pou5f1", "Pou3f1", "Utf1", "Hoxb1", "Hoxb2", "Hoxa1", "Cer1", "Krt8", "Hoxa7", "Hoxb6", "Hoxc5", "Tbx4", "Asz1", "Dazl", "Noto" ,"Tbx3", "Tbx3",
"Mcm4", "Mcm7", "Top2a", "Ube2c", "Pcan", "Rps12", "Rpl12", "Rpl13", "Rps13","Rpl14", "Rps14", "Xlr3b", "Aurkc", "Gstm6", "Fgf23", "Pnp2", "T", "Dnmt3b", "Dnmt3a", "Pim2", "Dab2", "Psme1", "Psme2", "Gata5", "Klf3", "Dppa3", "Cdkn2a")

init = function()
{
	mod = list()
	options(gmax.data.size=1e+9)
	options(scipen=999)

	mod$cov_stat = readRDS(file=in_path("data/output/multi_type_cov.RDS"))
	rownames(mod$cov_stat) = mod$cov_stat$nm
	mod$cov_peak_chr = readRDS(file=in_path("data/output/cov_peak_chr.RDS"))
	mod$cov_chr = readRDS(file=in_path("data/output/cov_chr.RDS"))
	mod$npeaks = readRDS(file=in_path("data/output/npeaks.RDS"))
	mod$cgweak_intervs = readRDS(file=in_path("data/output/cgweak_intervs.RDS"))

	mod$cov_states = c("marginal","rnp","exe_meso","ac_prog", "ac","haem_prog", "lnm","pmeso","card_meso","se","gut")
	mod$tss = gintervals.load("intervs.global.tss")
	mod$exon = gintervals.load("intervs.global.exon")
	mod = gen_cgi(mod)
	mod = gen_ctcf(mod)
	mod = init_exp(mod)

	mod$r_chr_cov_dnmt = mod$cov_chr[["DNMT_TKO_KO"]]/(mod$cov_chr[["WT_sm"]]+mod$cov_chr[["WT_lg"]])
	mod$r_chr_cov_dnmtd = mod$cov_chr[["DNMT_DKO_KO"]]/(mod$cov_chr[["WT_sm"]]+mod$cov_chr[["WT_lg"]])
	mod$r_chr_cov_dnmt_both = (mod$cov_chr[["DNMT_TKO_KO"]]+mod$cov_chr[["DNMT_DKO_KO"]])/(mod$cov_chr[["WT_sm"]]+mod$cov_chr[["WT_lg"]])
	mod$r_chr_cov_tet = mod$cov_chr[["TET_TKO_KO"]]/(mod$cov_chr[["WT_sm"]]+mod$cov_chr[["WT_lg"]])
	mod$r_chr_cov_dnmt = mod$r_chr_cov_dnmt/apply(mod$r_chr_cov_dnmt, 1, median, na.rm=T)
	mod$r_chr_cov_dnmtd = mod$r_chr_cov_dnmtd/apply(mod$r_chr_cov_dnmtd, 1, median, na.rm=T)
	mod$r_chr_cov_dnmt_both = mod$r_chr_cov_dnmt_both/apply(mod$r_chr_cov_dnmt_both, 1, median, na.rm=T)
	mod$r_chr_cov_tet = mod$r_chr_cov_tet/apply(mod$r_chr_cov_tet, 1, median, na.rm=T)

	mod$tet_tot_karyo = rowSums(mod$cov_chr[["TET_TKO_KO"]],na.rm=T)/rowSums(mod$cov_chr[["TET_TKO_KO"]]/mod$r_chr_cov_tet, na.rm=T)
	mod$dnmt_tot_karyo = rowSums(mod$cov_chr[["DNMT_TKO_KO"]],na.rm=T)/rowSums(mod$cov_chr[["DNMT_TKO_KO"]]/mod$r_chr_cov_dnmt, na.rm=T)
	mod$dnmtd_tot_karyo = rowSums(mod$cov_chr[["DNMT_DKO_KO"]],na.rm=T)/rowSums(mod$cov_chr[["DNMT_DKO_KO"]]/mod$r_chr_cov_dnmtd, na.rm=T)
	
	cov_norm = mod$cov_stat
	cov_norm$TET_TKO_KO.tot = floor(mod$cov_stat$TET_TKO_KO.tot / mod$tet_tot_karyo)
	cov_norm$TET_TKO_KO.tot_peak = floor(mod$cov_stat$TET_TKO_KO.tot_peak / mod$tet_tot_karyo)
	cov_norm$TET_TKO_KO.tot_rp = floor(mod$cov_stat$TET_TKO_KO.tot_rp / mod$tet_tot_karyo)

	cov_norm$DNMT_TKO_KO.tot = floor(mod$cov_stat$DNMT_TKO_KO.tot / mod$dnmt_tot_karyo)
	cov_norm$DNMT_TKO_KO.tot_peak = floor(mod$cov_stat$DNMT_TKO_KO.tot_peak / mod$dnmt_tot_karyo)
	cov_norm$DNMT_TKO_KO.tot_rp = floor(mod$cov_stat$DNMT_TKO_KO.tot_rp / mod$dnmt_tot_karyo)

	cov_norm$DNMT_DKO_KO.tot = floor(mod$cov_stat$DNMT_DKO_KO.tot / mod$dnmtd_tot_karyo)
	cov_norm$DNMT_DKO_KO.tot_peak = floor(mod$cov_stat$DNMT_DKO_KO.tot_peak / mod$dnmtd_tot_karyo)
	cov_norm$DNMT_DKO_KO.tot_rp = floor(mod$cov_stat$DNMT_DKO_KO.tot_rp / mod$dnmtd_tot_karyo)

	mod$cov_stat_kary_norm = cov_norm

	mod = update_downsamp_ratios(mod)

	# message("building ehs modes")
	if(!file.exists("data/output/epimod_ehs.RDS")) {
		mod = build_epi_hotspots_classes(mod)
		saveRDS(mod$ehs, out_path("data/output/epimod_ehs.RDS"))
	} else {
		mod$ehs = readRDS(in_path("data/output/epimod_ehs.RDS"))
	}

	mod$type_legc = readRDS(file=in_path("data/output/spatwt_type_legc.RDS"))

	type_key = read.csv("data/ref/cell_type_key.csv")
	mod$type_code = type_key$shorthand
	names(mod$type_code) = type_key$cell_type
	mod$code_type = type_key$cell_type
	names(mod$code_type) = type_key$shorthand
	mod$code_color = type_key$color
	names(mod$code_color) = type_key$shorthand

	if(file.exists("data/output/peak_stat.RDS")) {
		mod$peak_stat = readRDS(file=in_path("data/output/peak_stat.RDS"))
	} else {
		peak_stat = gen_peaks_z_mat(mod, mod$npeaks)
		mod$peak_stat = peak_stat
		saveRDS(peak_stat, file=out_path("data/output/peak_stat.RDS"))
	}
	return(mod)
}

init_exp = function(mod)
{
	col_fn = "data/ref/cell_type_colors.tsv"
	coltab = read.csv(col_fn, stringsAsFactors=F, sep='\t')
	type_col = coltab$color
	names(type_col) = coltab$cell_type

	mod$type_col = type_col

	if(!file.exists("data/output/spatwt_mc_legc.RDS")) {
		mc_fn = paste0(spatwt_dir, "/spatwt.metacells.h5ad")
		mcdata = read_h5ad(mc_fn)

		mc_egc = t(as.matrix(mcdata$X))
		mc_legc = log2(1e-5+mc_egc)
		mod$mc_legc =  mc_legc
		mod$mc_type = as.character(mcdata$obs$cell_type)
		names(mod$mc_type) = mcdata$obs_names

		type_legc = t(tgs_matrix_tapply(mod$mc_legc, mod$mc_type, mean))
		mod$type_legc = type_legc
		saveRDS(type_legc, file=out_path("data/output/spatwt_type_legc.RDS"))
		saveRDS(mc_legc, file=out_path("data/output/spatwt_mc_legc.RDS"))
	} else {
		mod$type_legc =  readRDS(file=in_path("data/output/spatwt_type_legc.RDS"))
		mod$mc_legc = readRDS(file=in_path("data/output/spatwt_mc_legc.RDS"))
	}

	return(mod)
}

gen_cgi = function(mod)
{

	if(!file.exists("data/output/cgi.RDS")) {	
		gvtrack.create("cg500", "seq.CG_new", "sum")
		gvtrack.iterator("cg500", sshift=-240, eshift=240)
		gvtrack.create("gc500", "seq.G_or_C", "sum")
		gvtrack.iterator("gc500", sshift=-240, eshift=240)
		cgi = gscreen("gc500>250 & (4*500*cg500/(gc500*gc500))>0.6", iterator=20)
		mod$cgi = cgi[(cgi$end-cgi$start) > 200,]
		saveRDS(cgi, file=out_path("data/output/cgi.RDS"))
	} else {
		mod$cgi = readRDS(in_path("data/output/cgi.RDS"))
	}
	return(mod)
}

gen_ctcf = function(mod)
{

	if(!file.exists("data/output/ctcf_peak.RDS")) {	
		ctcf = gscreen("chipseq.boyan.ES_CTCF_2reps>32")

		ctcf = ctcf[ctcf$chrom != "chrM" & ctcf$chrom != "chrY",]

		ctcf$start = ctcf$start - 100
		ctcf$end = ctcf$end + 100

		gvtrack.create("e_ctcf_f_max", "motifs.CTCF_forward1", "max")
		gvtrack.create("e_ctcf_r_max", "motifs.CTCF_reverse1", "max")

		ctcf_v = gextract(c("e_ctcf_f_max", "e_ctcf_r_max", "chipseq.boyan.ES_CTCF_2reps"), iterator=1, intervals=ctcf, colnames=c("e_f", "e_r", "chip"))
		max_chip = tapply(ctcf_v$chip, ctcf_v$intervalID, max)
		max_e_f = tapply(ctcf_v$e_f, ctcf_v$intervalID, max)
		max_e_r = tapply(ctcf_v$e_r, ctcf_v$intervalID, max)
		max_e = pmax(max_e_f, max_e_r)
		ctcf_peak_f = ctcf_v[ctcf_v$e_f==max_e[ctcf_v$intervalID],]
		ctcf_peak_r = ctcf_v[ctcf_v$e_r==max_e[ctcf_v$intervalID],]
		ctcf_peak_f$e_CTCF = ctcf_peak_f$e_f
		ctcf_peak_f$strand = 1
		ctcf_peak_r$e_CTCF = ctcf_peak_r$e_r
		ctcf_peak_r$strand = -1
		ctcf_peak = rbind(ctcf_peak_f, ctcf_peak_r)

		ctcf_peak = ctcf_peak[!duplicated(ctcf_peak$intervalID),]
		ctcf_peak$max_chip = max_chip[ctcf_peak$intervalID]

		saveRDS(ctcf_peak, file=out_path("data/output/ctcf_peak.RDS"))
		mod$ctcf_peak = ctcf_peak
	} else {
		mod$ctcf_peak = readRDS(in_path("data/output/ctcf_peak.RDS"))
	}
	return(mod)
}

gen_cgweak_interv = function()
{
    gvtrack.create("GC500", "seq.G_or_C", "sum")
	cgweak = gscreen("seq.CG_500_mean_new>0.022")

	gvtrack.create("e85_meth", "wgPBAT.E8_5.avg", "avg")
	gvtrack.iterator("e85_meth", eshift = -140, sshift=140)
	gvtrack.create("d_tss", "intervs.global.tss", "distance")
	cgweak_cg = gextract(c("seq.CG_500_mean_new", "GC500", "d_tss","e85_meth"), iterator=cgweak, 
						intervals=cgweak,colnames=c("CG","GC","d_tss","e85_meth"))
	cgweak_intervs = cgweak_cg
	cgweak_intervs = cgweak_intervs[cgweak_intervs$chrom != "chrY" & cgweak_intervs$chrom != "chrM",]
	saveRDS(cgweak_intervs, file=out_path("data/output/cgweak_intervs.RDS"))

}

build_epi_hotspots_classes = function(mod)
{
#merge intervals: peaks, high CG, CTCF
	mat_acp = build_peak_mat(mod, mod$npeaks, "ac_prog")
	mat_rnp = build_peak_mat(mod, mod$npeaks, "rnp")
	wt_rnp = rowSums(mat_rnp[,4:8])
	wt_acp = rowSums(mat_acp[,4:8])
	wt_rnp_p = wt_rnp/sum(wt_rnp)
	wt_acp_p = wt_acp/sum(wt_acp)
	mod$npeaks$rnp_p = log2(1e-6+wt_rnp_p)
	mod$npeaks$acp_p = log2(1e-6+wt_acp_p)

	tss = gintervals.load("intervs.global.tss")
	tss = tss[tss$chrom != "chrM" & tss$chrom != "chrY",]
    
	u1 = gintervals.union(intervals1=mod$npeaks, intervals2=mod$ctcf_peak)
	u2 = gintervals.union(intervals1=u1, intervals2 = tss)
	ehs = gintervals.union(intervals1=u2, intervals2 = mod$cgweak_intervs)

	ehs = ehs[ehs$chrom != "chrM" & ehs$chrom != "chrY",]

	gvtrack.create("cg_max", "seq.CG_500_mean_new", "max")
	e_cg = gextract(c("cg_max"), intervals=ehs, iterator=ehs, colnames=c("cg_max"))
	ehs$cg_max = e_cg$cg_max

	ehs_peak = gintervals.neighbors1(ehs, mod$npeaks, na.if.notfound=T)
	ehs_tss = gintervals.neighbors1(ehs, tss, na.if.notfound=T)
	ehs_cgweak = gintervals.neighbors1(ehs, mod$cgweak_intervs, na.if.notfound=T)
	ehs_ctcf = gintervals.neighbors1(ehs, mod$ctcf_peak, na.if.notfound=T)

	ehs$atac_d = ehs_peak$dist
	ehs$atac_rnp_p = ehs_peak$rnp_p
	ehs$atac_acp_p = ehs_peak$acp_p

	ehs$ctcf_d = ehs_ctcf$dist
	ehs$ctcf_chip = ehs_ctcf$chip
	ehs$ctcf_e = ehs_ctcf$e_CTCF
	ehs$ctcf_strand = ehs_ctcf$strand

	ehs$tss_d = ehs_tss$dist
	ehs$tss_strand = ehs_tss$strand

	ehs$cg_d = ehs_cgweak$dist
	ehs$cg_CG = ehs_cgweak$CG
	ehs$cg_GC = ehs_cgweak$GC

#partition all hotspts:
	ehs$mode = "none";
	ehs$center = floor((ehs$start+ehs$end)/2)
	ehs$strand = ehs$tss_strand

#CTCF non TSS
	f_ins_sctcf = ehs$ctcf_d == 0 & ehs$ctcf_chip >= 50 & abs(ehs$tss_d) >= 1000
	ehs$mode[f_ins_sctcf] = "ins_sctcf"
	f_ins_wctcf = ehs$ctcf_d == 0 & ehs$ctcf_chip < 50 & abs(ehs$tss_d) >= 1000
	ehs$mode[f_ins_wctcf] = "ins_wctcf"
	f_ins_ctcf = f_ins_sctcf | f_ins_wctcf
	ehs$strand[f_ins_ctcf] = ehs$ctcf_strand[f_ins_ctcf]
	ehs$center[f_ins_ctcf] = ehs_ctcf$start1[f_ins_ctcf]

#TSS w CTCF (center on CTCF)
	f_tss_ctcf = ehs$ctcf_d == 0 & abs(ehs$tss_d) < 1000
	ehs$mode[f_tss_ctcf] = "tss_ctcf"
	ehs$strand[f_tss_ctcf] = ehs$tss_strand[f_tss_ctcf]
	ehs$center[f_tss_ctcf] = ehs_ctcf$start1[f_tss_ctcf]

#TSS w/o CTCF
	f_tss = abs(ehs$tss_d) == 0 & ehs$mod == "none"
	ehs$mode[f_tss & ehs$cg_d == 0 & ehs$cg_max>0.04] = "tss_hcg"
	ehs$mode[f_tss & ehs$cg_d == 0 & ehs$cg_max<0.04] = "tss_icg"
	ehs$mode[f_tss & ehs$cg_d != 0] = "tss_lcg"
	ehs$strand[f_tss] = ehs$tss_strand[f_tss]
	ehs$center[f_tss] = ehs_tss$start1[f_tss]

#atac high CG non TSS
	atac_p = (ehs$atac_rnp_p+ehs$atac_acp_p)/2
	ehs$mode[abs(ehs$cg_d) == 0 & ehs$cg_max > 0.04 & abs(ehs$atac_d) == 0 & atac_p >= -16.5 & ehs$mod == "none"] = "hcg_senh"
	ehs$mode[abs(ehs$cg_d) == 0 & ehs$cg_max > 0.04 & abs(ehs$atac_d) == 0 & atac_p < -16.5 & ehs$mod == "none"] = "hcg_wenh"
#atac mid CG non TSS
	ehs$mode[abs(ehs$cg_d) == 0 & ehs$cg_max < 0.04 & abs(ehs$atac_d) == 0 & atac_p >= -16.5 & ehs$mod == "none"] = "icg_senh"
	ehs$mode[abs(ehs$cg_d) == 0 & ehs$cg_max < 0.04 & abs(ehs$atac_d) == 0 & atac_p < -16.5 & ehs$mod == "none"] = "icg_wenh"
#atac low CG non TSS
	ehs$mode[ehs$atac_d == 0 & atac_p >= -16.5 & ehs$mod == "none"] = "lcg_senh"
	ehs$mode[ehs$atac_d == 0 & atac_p >= -18.5 & ehs$mod == "none"] = "lcg_menh"
	ehs$mode[ehs$atac_d == 0 & atac_p < -18.5 & ehs$mod == "none"] = "lcg_wenh"

#atac high CG non TSS
	ehs$mode[abs(ehs$cg_d) == 0 & ehs$cg_max > 0.04 & ehs$mod == "none"] = "hcg_void"
#atac mid CG non TSS
	ehs$mode[abs(ehs$cg_d) == 0 & ehs$cg_max < 0.04 & ehs$mod == "none"] = "icg_void"
#atac low CG non TSS
    ehs$mode[abs(ehs$cg_d) != 0 & ehs$cg_max < 0.04 & ehs$mod == "none"] = "lcg_void"

	mod$ehs = ehs
	return(mod)
}

gen_track_vext = function(tname = "marginal", w_ext = 140)
{
	message("set vt for ", tname)
	nm_dnmt = sprintf("mut_multi2.DNMT_TKO_KO.%s", tname)
	nm_dnmtd = sprintf("mut_multi2.DNMT_DKO_KO.%s", tname)
	nm_tet = sprintf("mut_multi2.TET_TKO_KO.%s", tname)
	nm_wtl = sprintf("mut_multi2.WT_lg.%s", tname)
	nm_wts = sprintf("mut_multi2.WT_sm.%s", tname)
	nm_wtht = sprintf("mut_multi2.DNMT_TKO_host.%s", tname)
	nm_wthd = sprintf("mut_multi2.DNMT_DKO_host.%s", tname)
	nm_wthtet = sprintf("mut_multi2.TET_TKO_host.%s", tname)
	gvtrack.create("dnmt_ext", nm_dnmt, "sum")
	gvtrack.create("dnmtd_ext", nm_dnmtd, "sum")
	gvtrack.create("tet_ext", nm_tet, "sum")
	gvtrack.create("wtl_ext", nm_wtl, "sum")
	gvtrack.create("wts_ext", nm_wts, "sum")
	gvtrack.create("wtht_ext", nm_wtht, "sum")
	gvtrack.create("wthd_ext", nm_wthd, "sum")
	gvtrack.create("wthtet_ext", nm_wthtet, "sum")
	gvtrack.iterator("dnmt_ext", sshift = -w_ext, eshift = w_ext)
	gvtrack.iterator("dnmtd_ext", sshift = -w_ext, eshift = w_ext)
	gvtrack.iterator("tet_ext", sshift = -w_ext, eshift = w_ext)
	gvtrack.iterator("wtl_ext", sshift = -w_ext, eshift = w_ext)
	gvtrack.iterator("wts_ext", sshift = -w_ext, eshift = w_ext)
	gvtrack.iterator("wtht_ext", sshift = -w_ext, eshift = w_ext)
	gvtrack.iterator("wthd_ext", sshift = -w_ext, eshift = w_ext)
	gvtrack.iterator("wthtet_ext", sshift = -w_ext, eshift = w_ext)
}

collect_cov_r = function()
{
	npeaks = readRDS(in_path("data/output/npeaks.RDS"))
	tss = gintervals.load("intervs.global.tss")
	tss = tss[tss$chrom != "chrM" & tss$chrom != "chrY",]

	all_nms = gtrack.ls("mut_multi2.WT_sm.")
	all_nms = union(all_nms, gtrack.ls("mut_multi2.WT_lg."))
	all_nms = union(all_nms, gtrack.ls("mut_multi2.DNMT_TKO_KO."))
	all_nms = sub("mut_multi2.WT_sm.","", all_nms)
	all_nms = sub("mut_multi2.WT_lg.","", all_nms)
	all_nms = sub("mut_multi2.DNMT_TKO_KO.","", all_nms)
	all_nms = all_nms[!grepl("_normed", all_nms) & !grepl("_smoothed", all_nms)]
	all_nms = unique(all_nms)

	f_rp = grepl("^Rp[l|s]", tss$geneSymbol)
	rp_tss = tss[f_rp,c("chrom","start", "end","strand")]
	rp_tss$start = rp_tss$start + ifelse(rp_tss$strand == 1, -800, -200)
	rp_tss$end = rp_tss$end + ifelse(rp_tss$strand == 1, 200, 800)

	cov_stats = list()
	cov_peak_chr = list()
	cov_chr = list()
	for(ds in c("TET_TKO_host", "DNMT_TKO_KO", "DNMT_DKO_KO", "TET_TKO_KO", "WT_lg", "WT_sm", "DNMT_TKO_host", "DNMT_DKO_host")) {
		cov_stat = NULL

		for(nm in all_nms) {
			message("doing track ", ds, " ", nm)
			track_nm = sprintf("mut_multi2.%s.%s", ds, nm)
			if(gtrack.exists(track_nm)) {
				tot = as.numeric(gsummary(track_nm)[["Sum"]])
				a = gextract(track_nm, intervals=npeaks, iterator=20, colnames=c("cov"))
				a_rp = gextract(track_nm, intervals=rp_tss, iterator=20, colnames=c("cov"))
				gvtrack.create("tcov", track_nm, "sum")
				a_chrom = gextract(c("tcov"), intervals=gintervals.all(), iterator=gintervals.all(), colnames=c("cov"))
				tot_peak = sum(a$cov, na.rm=T)
				tot_peak_chr = tapply(a$cov, a$chrom, sum, na.rm=T)
				tot_rp = sum(a_rp$cov, na.rm=T)
				if(is.null(cov_stat)) {
					cov_stat = data.frame("nm" =nm ,"tot"=tot, "tot_peak"=tot_peak, "tot_rp" = tot_rp)
					cov_peak_chr_st = c(nm, tot_peak_chr)
					cov_chr_st = c(nm, a_chrom$cov)
				} else {
					cov_stat = rbind(cov_stat, c(nm, tot, tot_peak, tot_rp))		
					cov_peak_chr_st = rbind(cov_peak_chr_st, c(nm, tot_peak_chr))
					cov_chr_st = rbind(cov_chr_st, c(nm, a_chrom$cov))
				}
			} else {
				if(is.null(cov_stat)) {
					cov_stat = data.frame("nm"= nm, "tot" = NA, "tot_peak" = NA, "tot_rp" = NA)		
					cov_peak_chr_st = c(nm,rep(NA,22))
					cov_chr_st = c(nm,rep(NA,22))
				} else {
					cov_stat = rbind(cov_stat, c(nm, NA, NA, NA))		
					cov_peak_chr_st = rbind(cov_peak_chr_st, c(nm,rep(NA,22)))
					cov_chr_st = rbind(cov_chr_st, c(nm,rep(NA,22)))
				}
			}
		}
		cov_stat$tot = as.numeric(cov_stat$tot)
		cov_stat$tot_peak = as.numeric(cov_stat$tot_peak)
		colnames(cov_stat)[2] = sprintf("%s.tot", ds)
		colnames(cov_stat)[3] = sprintf("%s.tot_peak", ds)
		colnames(cov_stat)[4] = sprintf("%s.tot_rp", ds)
		cov_stats[[ds]] = cov_stat
		rownames(cov_peak_chr_st) = cov_peak_chr_st[,1]
		cov_peak_chr_st = cov_peak_chr_st[,-1]
		class(cov_peak_chr_st) = 'numeric'
		cov_peak_chr[[ds]] = cov_peak_chr_st

		rownames(cov_chr_st) = cov_chr_st[,1]
		cov_chr_st = cov_chr_st[,-1]
		colnames(cov_chr_st) = gintervals.all()$chrom
		class(cov_chr_st) = 'numeric'
		cov_chr[[ds]] = cov_chr_st
	}
	cov_stat = cov_stats[["DNMT_TKO_KO"]] %>% 
									left_join(cov_stats[["DNMT_DKO_KO"]]) %>%
									left_join(cov_stats[["TET_TKO_KO"]]) %>%
									left_join(cov_stats[["WT_lg"]]) %>%
									left_join(cov_stats[["WT_sm"]]) %>%
									left_join(cov_stats[["DNMT_TKO_host"]]) %>%
									left_join(cov_stats[["DNMT_DKO_host"]]) %>%
									left_join(cov_stats[["TET_TKO_host"]])

	cov_stat$WT_sm.tot_rp = as.numeric(cov_stat$WT_sm.tot_rp)
	cov_stat$WT_lg.tot_rp = as.numeric(cov_stat$WT_lg.tot_rp)
	cov_stat$TET_TKO_KO.tot_rp = as.numeric(cov_stat$TET_TKO_KO.tot_rp)
	cov_stat$DNMT_TKO_KO.tot_rp = as.numeric(cov_stat$DNMT_TKO_KO.tot_rp)
	cov_stat$DNMT_DKO_KO.tot_rp = as.numeric(cov_stat$DNMT_DKO_KO.tot_rp)
	cov_stat$DNMT_TKO_host.tot_rp = as.numeric(cov_stat$DNMT_TKO_host.tot_rp)
	cov_stat$DNMT_DKO_host.tot_rp = as.numeric(cov_stat$DNMT_DKO_host.tot_rp)
	cov_stat$TET_TKO_host.tot_rp = as.numeric(cov_stat$TET_TKO_host.tot_rp)

	saveRDS(cov_stat, file=out_path("data/output/multi_type_cov.RDS"))
	saveRDS(cov_peak_chr, file=out_path("data/output/cov_peak_chr.RDS"))
	saveRDS(cov_chr, file=out_path("data/output/cov_chr.RDS"))
	rownames(cov_stat) = cov_stat$nm

}

update_downsamp_ratios = function(mod) 
{
	tot_wt = rowSums(mod$cov_stat[,c("WT_lg.tot_rp","WT_sm.tot_rp", "DNMT_TKO_host.tot_rp","DNMT_DKO_host.tot_rp", "TET_TKO_host.tot_rp")], na.rm=T)

	tot_dnmt = mod$cov_stat[,"DNMT_TKO_KO.tot_rp"]
	tot_dnmtd = mod$cov_stat[,"DNMT_DKO_KO.tot_rp"]
	tot_dnmt_both = rowSums(mod$cov_stat[,c("DNMT_TKO_KO.tot_rp","DNMT_DKO_KO.tot_rp")], na.rm=T)
	tot_tet = mod$cov_stat[,"TET_TKO_KO.tot_rp"]

	min_ignore = pmin(pmin(tot_wt, tot_dnmt), tot_tet)
	min_compare = pmin(pmin(pmin(tot_wt, tot_dnmt), tot_dnmtd),tot_tet)
	min_pool = pmin(pmin(tot_wt, tot_dnmt_both),tot_tet)

	mod$r_to_dsamp[["ignore"]] = data.frame(wt = min_ignore/tot_wt, dnmt = min_ignore/tot_dnmt, tet = min_ignore/tot_tet)
	mod$r_to_dsamp[["compare"]] = data.frame(wt = min_compare/tot_wt, dnmt = min_compare/tot_dnmt, dnmtd = min_compare/tot_dnmtd, tet = min_compare/tot_tet)
	mod$r_to_dsamp[["pool"]] = data.frame(wt = min_pool/tot_wt, dnmt = min_pool/tot_dnmt_both, tet = min_pool/tot_tet)
	return(mod)
}

downsamp_karyo_balance = function(df, col_nm, karyo_r)
{
	bad_chr = intersect(names(which(karyo_r>1.06)), unique(as.character(df$chrom)))
	for(chr in bad_chr) {
		f = df$chrom == chr
		r = karyo_r[chr]
		mat = matrix(df[f, col_nm],ncol=1)
		message("correct karyo ", col_nm, " chr ", chr, " r = ",r, " n =", round(sum(mat)/r))
		mat_ds = downsample_matrix(mat,  target_n=round(sum(mat)/r), seed=19)
		df[f,col_nm] = mat_ds
	}
	return(df)
}

comp_capped_regional_atac = function(track_nm, q=0.95, win_local = 300, win_regional = 20000)
{
	gvtrack.create("cov_local", track_nm, "sum")
	local_ext = (win_local-20)/2
	local_n = win_local/20
	gvtrack.iterator("cov_local", sshift = -local_ext, eshift = local_ext)

	T_q = gquantiles("cov_local",q)

	cap_q = function(x) { x[is.na(x)] = 0; return(ifelse(x > T_q, T_q, x)/local_n) }
	gtrack.smooth(sprintf("%s_regcap", track_nm), 
						desc= "capped regional avg", 
						expr = "cap_q(cov_local)",
						winsize = win_regional,
						iterator=20)

}

comp_peaks_over_punc_region = function(track_nm, norm_track_nm, peak_ext, eps=3, T_q = 0.98)
{
	peak_l = peak_ext*2+20
	peak_n = peak_l/20
	gvtrack.create("cov_loc", track_nm, "sum")
	gvtrack.iterator("cov_loc", sshift=-peak_ext, eshift=peak_ext)	
	norm_expr = sprintf("(%s+cov_loc)/(%s+peak_n*%s)", eps, eps, norm_track_nm)
	T_rat = gquantiles(norm_expr, T_q, iterator=20)
	message("ratio thresh ", T_rat)
	peaks = gscreen(sprintf("%s > %s", norm_expr, T_rat), iterator=20)

	return(peaks)
}

join_norm_peaks = function(T_q = 0.97, eps=3, state = "marginal")
{
	base = "mut_multi2"
	norm_suf = "marginal_regcap"

	peaks_wtl = comp_peaks_over_punc_region(sprintf("%s.WT_lg.%s", base, state), 
												sprintf("%s.WT_lg.%s", base, norm_suf),
												peak_ext=140, eps=eps, T_q=T_q)
	peaks_wts = comp_peaks_over_punc_region(sprintf("%s.WT_sm.%s", base, state), 
												sprintf("%s.WT_sm.%s", base, norm_suf),
												peak_ext=140, eps=eps, T_q=T_q)
	WT_peaks = gintervals.union(peaks_wts, peaks_wtl)
	peaks_dnmt = comp_peaks_over_punc_region(sprintf("%s.DNMT_TKO_KO.%s", base, state), 
												sprintf("%s.DNMT_TKO_KO.%s", base, norm_suf),
												peak_ext=140, eps=eps, T_q=T_q)
	peaks_dnmtd = comp_peaks_over_punc_region(sprintf("%s.DNMT_DKO_KO.%s", base, state), 
												sprintf("%s.DNMT_DKO_KO.%s", base, norm_suf),
												peak_ext=140, eps=eps, T_q=T_q)
	peaks_tet = comp_peaks_over_punc_region(sprintf("%s.TET_TKO_KO.%s", base, state), 
												sprintf("%s.TET_TKO_KO.%s", base, norm_suf),
												peak_ext=140, eps=eps, T_q=T_q)
	mut_peaks = gintervals.union(gintervals.union(peaks_dnmt, peaks_dnmtd,), peaks_tet)

	peaks = gintervals.union(WT_peaks, mut_peaks)

	return(peaks)
}

build_peak_mat = function(mod, peaks, state, ignore_dko_dsamp = F)
{
	gen_track_vext(state, w_ext=140)

	bg_wtl = "mut_multi2.WT_lg.marginal_regcap"
	bg_wts = "mut_multi2.WT_sm.marginal_regcap"
	bg_dnmt = "mut_multi2.DNMT_TKO_KO.marginal_regcap"
	bg_tet = "mut_multi2.TET_TKO_KO.marginal_regcap"
	mat = gextract(c("wtl_ext", "wts_ext", "wthd_ext", "wtht_ext", "wthtet_ext", "dnmt_ext", "dnmtd_ext", "tet_ext", bg_wtl, bg_wts, bg_dnmt, bg_tet),
						intervals = peaks, iterator = peaks, 
						colnames = c("wtl_n", "wts_n", "wthd_n", "wtht_n", "wthtet_n", "dnmt_n", "dnmtd_n", "tet_n", "bg_wtl", "bg_wts", "bg_dnmt", "bg_tet"))

	mat[is.na(mat)] = 0

#correct karyotype
	mat = downsamp_karyo_balance(mat,"tet_n", mod$r_chr_cov_tet[state,])
	mat = downsamp_karyo_balance(mat,"dnmt_n", mod$r_chr_cov_dnmt[state,])
	mat = downsamp_karyo_balance(mat,"dnmtd_n", mod$r_chr_cov_dnmtd[state,])

	mat$wt_n = mat$wtl_n+mat$wts_n

#decide on ds coverage (these number already account for karyo correction)
	r_wt = mod$cov_stat_kary_norm[state,"WT_lg.tot_rp"] + 
								mod$cov_stat_kary_norm[state,"WT_sm.tot_rp"]
	r_dnmt = mod$cov_stat_kary_norm[state,"DNMT_TKO_KO.tot_rp"]
	r_dnmtd = mod$cov_stat_kary_norm[state,"DNMT_DKO_KO.tot_rp"]
	r_tet = mod$cov_stat_kary_norm[state,"TET_TKO_KO.tot_rp"]
	if(!ignore_dko_dsamp) {
		r_min = min(r_wt, r_dnmt, r_dnmtd, r_tet)
		r_dnmtd = r_min/r_dnmtd
	} else {
		r_min = min(r_wt, r_dnmt, r_tet)
		r_dnmtd = 1
	}
	r_wt = r_min/r_wt
	r_dnmt = r_min/r_dnmt
	r_tet = r_min/r_tet

#correct standartized coverage

	if(r_wt < 1) {
		to_ds = matrix(mat$wt_n,ncol=1)
		post_ds = downsample_matrix(to_ds,  target_n = floor(sum(to_ds)*r_wt), seed=19)
		mat$wt_ds = post_ds
	} else {
		mat$wt_ds = mat$wt_n
	}
	if(r_dnmt < 1) {
		to_ds = matrix(mat$dnmt_n,ncol=1)
		post_ds = downsample_matrix(to_ds,  target_n = floor(sum(to_ds)*r_dnmt), seed=19)
		mat$dnmt_ds = post_ds
	} else {
		mat$dnmt_ds = mat$dnmt_n
	}

	if(!ignore_dko_dsamp & r_dnmtd < 1) {
		to_ds = matrix(mat$dnmtd_n,ncol=1)
		post_ds = downsample_matrix(to_ds,  target_n = floor(sum(to_ds)*r_dnmtd), seed=19)
		mat$dnmtd_ds = post_ds
	} else {
		mat$dnmtd_ds = mat$dnmtd_n
	}
	if(r_tet < 1) {
		to_ds = matrix(mat$tet_n,ncol=1)
		post_ds = downsample_matrix(to_ds,  target_n = floor(sum(to_ds)*r_tet), seed=19)
		mat$tet_ds = post_ds
	} else {
		mat$tet_ds = mat$tet_n
	}
	mat$wt = mat$wt_n * r_wt
	mat$dnmt = mat$dnmt_n * r_dnmt
	mat$dnmtd = mat$dnmtd_n * r_dnmtd
	mat$tet= mat$tet_n * r_tet

	return(mat)	
}

gen_peaks_z_mat = function(mod, peaks, with_dko = T)
{
	rownames(peaks) = paste(peaks$chrom, peaks$start, sep=":")
	peak_p = list()
	for(state in mod$cov_states) {
		peak_p[[state]] = peaks
		message("extract data for ", state)
		mat = build_peak_mat(mod, peaks, state, ignore_dko_dsamp = !with_dko)
		z_dnmt = ifelse(mat$wt_ds + mat$dnmt_ds > 15, (mat$dnmt_ds - mat$wt_ds)/sqrt((mat$wt_ds + mat$dnmt_ds)/2),NA)
		r_dnmt = ifelse(mat$wt_ds + mat$dnmt_ds > 15, log2((mat$dnmt_ds+1)/(mat$wt_ds+1)),NA)
		peak_p[[state]]$z_dnmt = z_dnmt
		peak_p[[state]]$r_dnmt = r_dnmt

		z_tet = ifelse(mat$wt_ds + mat$tet_ds > 15, (mat$tet_ds - mat$wt_ds)/sqrt((mat$wt_ds + mat$tet_ds)/2),NA)
		r_tet = ifelse(mat$wt_ds + mat$tet_ds > 15, log2((mat$tet_ds+1)/(mat$wt_ds+1)),NA)
		peak_p[[state]]$z_tet = z_tet
		peak_p[[state]]$r_tet = r_tet
		z_dnte = ifelse(mat$dnmt_ds + mat$tet_ds > 15, (mat$dnmt_ds - mat$tet_ds)/sqrt((mat$dnmt_ds + mat$tet_ds)/2),NA)
		r_dnte = ifelse(mat$dnmt_ds + mat$tet_ds > 15, log2((mat$dnmt_ds+1)/(mat$tet_ds+1)),NA)
		peak_p[[state]]$z_dnte = z_dnte
		peak_p[[state]]$r_dnte = r_dnte
		peak_p[[state]]$dnmt_ds = mat$dnmt_ds
		peak_p[[state]]$tet_ds = mat$tet_ds
		peak_p[[state]]$wt_ds = mat$wt_ds
		peak_p[[state]]$dnmt_n = mat$dnmt_n
		peak_p[[state]]$tet_n = mat$tet_n
		peak_p[[state]]$wts_n = mat$wts_n
		peak_p[[state]]$wtl_n = mat$wtl_n

		if(with_dko) {
			z_dnmtd = ifelse(mat$wt_ds + mat$dnmtd_ds > 15, (mat$dnmtd_ds - mat$wt_ds)/sqrt((mat$wt_ds + mat$dnmtd_ds)/2),NA)
			r_dnmtd = ifelse(mat$wt_ds + mat$dnmtd_ds > 15, log2((mat$dnmtd_ds+1)/(mat$wt_ds+1)),NA)
			peak_p[[state]]$z_dnmtd = z_dnmtd
			peak_p[[state]]$r_dnmtd = r_dnmtd
			peak_p[[state]]$dnmtd_ds = mat$dnmtd_ds
		}
		peak_p[[state]]$dnmtd_n = mat$dnmtd_n
	}
	return(peak_p)
}

# ---- used by ExtData_ATAC.ipynb ----

add_mean_to_previous <- function(x) {
    n <- length(x)
    result <- numeric(n)

    for (i in 1:n) {
        if (i == 1) {
            result[i] <- x[i] / 2
        } else {
            diff <- (x[i] - x[i-1])
            result[i] <- x[i-1] + diff/2
        }
    }
  
    return(result)
}
