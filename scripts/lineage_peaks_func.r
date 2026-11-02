library(prego)
library(misha)
library(tibble)
library(ggplot2)
library(purrr)

gsetroot(mm10_trackdb)

cell_type_key = read.csv("./data/ref/cell_type_key.csv", row.names = 1)

gen_peaks_lin = function(mod)
{
	if(file.exists("data/output/lin_peaks.RDS")) {
		mod$lin_peaks = readRDS(in_path("data/output/lin_peaks.RDS"))
		return(mod)
	}

	# If mod$npeaks is missing columns CG, GC, or d_tss, add them
	if (!all(c("CG", "GC", "d_tss") %in% colnames(mod$npeaks))) {
		gvtrack.create("GC500", "seq.G_or_C", "sum")
		gvtrack.create("CG_max", "seq.CG_500_mean_new", "max")
		gvtrack.create("d_tss", "intervs.global.tss", "distance")
		ge = gextract(c("CG_max", "GC500/500", "d_tss"), intervals = mod$npeaks, iterator = mod$npeaks, colnames = c("CG", "GC", "d_tss"))
		mod$npeaks$CG = ge$CG
		mod$npeaks$GC = ge$GC
		mod$npeaks$d_tss = ge$d_tss
	}

	pms = list()
	states = c("rnp", "ac_prog", "ac", "pmeso", "card_meso", "lnm", "ery1","se")
	for(st in states) {
 	     pms[[st]] = build_peak_mat(mod, mod$npeaks, st)
	}

	all = lapply(pms, function(x) return(cbind(x[,4]+x[,5],x[,c(9,10,11)])))
	mall = do.call('cbind', all)
	colnames(mall)[c(1,5,9,13,17,21,25,29)] = paste(states, "wt", sep=".")

	rnp = mall[,1:4]
	acp = mall[,5:8]+mall[,9:12] #acp+ac
	meso = mall[,13:16]+mall[,17:20] #pmeso+cardmeso+lnm
	lnm = mall[,21:24] #pmeso+cardmeso+lnm
	ery = mall[,25:28] #ery1
	se = mall[,29:32] #se

	rnp_ds = downsample_matrix(as.matrix(rnp),target_n=min(colSums(rnp[,-3])))
	acp_ds = downsample_matrix(as.matrix(acp),target_n=min(colSums(acp[,-3])))
	meso_ds = downsample_matrix(as.matrix(meso), target_n=min(colSums(meso[,-3])))
	lnm_ds = downsample_matrix(as.matrix(lnm), target_n=min(colSums(lnm[,-3])))
	ery_ds = downsample_matrix(as.matrix(ery), target_n=min(colSums(ery[,-3])))
	se_ds = downsample_matrix(as.matrix(se), target_n=min(colSums(se[,-3])))

	mod$lin_peaks = list(rnp_ds = rnp_ds, acp_ds = acp_ds, meso_ds = meso_ds, lnm_ds = lnm_ds, ery_ds = ery_ds, se_ds = se_ds, m_all = mall)
	saveRDS(mod$lin_peaks, out_path("data/output/lin_peaks.RDS"))
	return(mod)
}

plot_trends_cgbinned <- function(mod){

	# Bin by CG, GC content
	cgbin = as.numeric(cut(mod$npeaks$CG,c(seq(0,0.05,0.01),1), include.lowest = T))
	gcbin = as.numeric(cut(mod$npeaks$GC,c(0,seq(0.4,0.6,0.05), include.lowest = T)))

	f_tss = abs(mod$npeaks$d_tss)<1000 | abs(mod$npeaks$d_tss)>1e+6
	f_hgc = mod$npeaks$GC > 0.55
	f_lcg = mod$npeaks$CG < 0.015

	# Compute deltas
	atac1 =  mod$lin_peaks$acp_ds
	atac2 = mod$lin_peaks$rnp_ds

	e1 = t(t(atac1)/colSums(atac1))
	e2 = t(t(atac2)/colSums(atac2))
	le1 = log2(1e-6+e1)
	le2 = log2(1e-6+e2)

	dlt1 = le1[,2]-le1[,4]
	dlt2 = le2[,2]-le2[,4]

	# ACP: elements per bin barplots, delta boxplots for DNMT, TET
	plot_df = le1 %>% as.data.frame() %>% mutate(WT_bin = ifelse(ac_prog.wt < -16.5, "weak", "strong"), 
                                                 mean = apply(le1[,c(1,2,4)], 1, mean),
                                                 mean_bin = ifelse(mean < -16.5, "weak", "strong"), 
                                                 CG_bin = as.character(cgbin))

	plot_df = plot_df[!f_tss,] %>% filter(!is.na(CG_bin)) %>% mutate(CG_bin = factor(CG_bin))
	levels(plot_df$CG_bin) = c(paste0("<0.0", 1:5), ">0.05")
	plot_df = plot_df %>% filter(!is.na(CG_bin))
	plot_df = plot_df %>% filter(!is.na(CG_bin))

	plot_df %>% group_by(mean_bin) %>% mutate(N = n()) %>% ungroup() %>% group_by(mean_bin, CG_bin) %>% summarize(f = n()/N) %>% ungroup() %>%
		ggplot(aes(CG_bin, f, fill = mean_bin)) + geom_bar(width = 0.8, stat="identity", position="dodge") +
			scale_fill_manual(values = c("weak" = "#5DD6FF", "strong" = "deepskyblue3")) + theme_minimal() +
			theme(legend.position = "None", strip.text.x = element_text(size = 15)) + xlab("%CG") + ylab("Fraction of Loci") + facet_wrap(~mean_bin)
	show_fig("./figs/Fig3/ACP_n_barplots.pdf", ggplot2::last_plot(), width = 6, height = 2)

	ggplot(plot_df, aes(CG_bin, ac_prog.dnmt_n-ac_prog.wt, fill = mean_bin)) + geom_boxplot(width = 0.5, outlier.size = 0.5) + ylim(-1.5,1.5) +
		scale_fill_manual(values = c("weak" = "#5DD6FF", "strong" = "deepskyblue3")) + theme_minimal() +
		xlab("%CG") + ylab("ACP") + theme(legend.position = "None", strip.text.x = element_text(size = 15)) + facet_wrap(~mean_bin) + geom_hline(yintercept = 0)
	show_fig("./figs/Fig3/DNMT_ACP_boxplots.pdf", ggplot2::last_plot(), width = 6, height = 2)
	ggplot(plot_df, aes(CG_bin, ac_prog.tet_n-ac_prog.wt, fill = mean_bin)) + geom_boxplot(width = 0.5, outlier.size = 0.5) + ylim(-1.5,1.5) +
		scale_fill_manual(values = c("weak" = "#5DD6FF", "strong" = "deepskyblue3")) + theme_minimal() +
		xlab("%CG") + ylab("ACP") + theme(legend.position = "None", strip.text.x = element_text(size = 15)) + facet_wrap(~mean_bin) + geom_hline(yintercept = 0)
	show_fig("./figs/Fig3/TET_ACP_boxplots.pdf", ggplot2::last_plot(), width = 6, height = 2)
    
    saveRDS(plot_df, out_path(paste0(embdnmt_dir, "/data/ref/Fig3d_ACP_df.RDS")))

	# RNP: elements per bin barplots, delta boxplots for DNMT, TET
	plot_df = le2 %>% as.data.frame() %>% mutate(WT_bin = ifelse(rnp.wt < -16.5, "weak", "strong"), 
												mean = apply(le2[,c(1,2,4)], 1, mean),
												mean_bin = ifelse(mean < -16.5, "weak", "strong"), 
												CG_bin = as.character(cgbin))
	plot_df = plot_df[!f_tss,] %>% filter(!is.na(CG_bin)) %>% mutate(CG_bin = factor(CG_bin))
	levels(plot_df$CG_bin) = c(paste0("<0.0", 1:5), ">0.05")
	plot_df = plot_df %>% filter(!is.na(CG_bin))

	plot_df %>% group_by(mean_bin) %>% mutate(N = n()) %>% ungroup() %>% group_by(mean_bin, CG_bin) %>% summarize(f = n()/N) %>% ungroup() %>%
		ggplot(aes(CG_bin, f, fill = mean_bin)) + geom_bar(width = 0.8, stat="identity", position="dodge") +
			scale_fill_manual(values = c("weak" = "#5DD6FF", "strong" = "deepskyblue3")) + theme_minimal() +
			theme(legend.position = "None", strip.text.x = element_text(size = 15)) + xlab("%CG") + ylab("Fraction of Loci") + facet_wrap(~mean_bin)
	show_fig("./figs/Fig3/RNP_n_barplots.pdf", ggplot2::last_plot(), width = 6, height = 2)

	ggplot(plot_df, aes(CG_bin, rnp.dnmt_n-rnp.wt, fill = mean_bin)) + geom_boxplot(width = 0.5, outlier.size = 0.5) + ylim(-1.5,1.5) +
		scale_fill_manual(values = c("weak" = "#5DD6FF", "strong" = "deepskyblue3")) + theme_minimal() +
		xlab("%CG") + ylab("RNP") + theme(legend.position = "None", strip.text.x = element_text(size = 15)) + facet_wrap(~mean_bin) + geom_hline(yintercept = 0)
	show_fig("./figs/Fig3/DNMT_RNP_boxplots.pdf", ggplot2::last_plot(), width = 6, height = 2)
	ggplot(plot_df, aes(CG_bin, rnp.tet_n-rnp.wt, fill = mean_bin)) + geom_boxplot(width = 0.5, outlier.size = 0.5) + ylim(-1.5,1.5) +
		scale_fill_manual(values = c("weak" = "#5DD6FF", "strong" = "deepskyblue3")) + theme_minimal() +
		xlab("%CG") + ylab("RNP") + theme(legend.position = "None", strip.text.x = element_text(size = 15)) + facet_wrap(~mean_bin) + geom_hline(yintercept = 0)
	show_fig("./figs/Fig3/TET_RNP_boxplots.pdf", ggplot2::last_plot(), width = 6, height = 2)
    
    saveRDS(plot_df, out_path(paste0(embdnmt_dir, "/data/ref/Fig3d_RNP_df.RDS")))

}

plot_scatters_cgbinned <- function(mod, suffix, nm1 = "rnp", nm2 = "acp", 
                                   strong_modes = c("hcg_senh", "icg_senh", "lcg_senh"), 
                                   weak_modes = c("hcg_wenh", "icg_wenh", "lcg_wenh"), 
                                   background_color = "grey70", weak_color = "#f7934c", strong_color = "#cc5803", 
                                   binned_by = "CG", save_RDS = F, add_labels = F){
    
    lp = readRDS("data/output/lin_peaks_ehs.RDS")
	bin = as.numeric(cut(mod$ehs[,binned_by] ,c(seq(0,0.05,0.01),1)))
    
    mode_vec = mod$ehs$mode
    chosen_modes = c(weak_modes, strong_modes)
    filter = mode_vec %in% chosen_modes
    
    geneSymbol_vec = gintervals.neighbors(mod$ehs, mod$tss, mindist = 0, maxdist = 0, na.if.notfound = T) %>% pull(geneSymbol)
    germline_genes = read.csv("./data/download/weber_genes.csv") %>% filter(class == "upregulated_germline") %>% pull(gene)
    
    atac1 =  lp$rnp_ds # RNP
	e1 = t(t(atac1)/colSums(atac1))
	le1 = log2(1e-6+e1)

	show_fig(sprintf("figs/Fig3/binned_scatter_%s_%s%s.png", binned_by, nm1, suffix), height = 1800, width = 5400, res = 300, draw = {
	layout(matrix(1:12, nrow=2))
	for(i in 1:6) {
        par(mar=c(3, 3, 2, 1))
		f = bin==i & filter
        f_up = (le1[f,2]-le1[f,1]) > 1
        f_down = (le1[f,2]-le1[f,1]) < -1
        f_weak = mode_vec[f] %in% weak_modes
        f_strong = mode_vec[f] %in% strong_modes
        
        f_germline = geneSymbol_vec[f] %in% germline_genes
        f_chrx =
        
        
        
		up = round(100*sum(f_up,na.rm=T)/sum(f,na.rm=T), 1)
		down = round(100*sum(f_down,na.rm=T)/sum(f,na.rm=T), 1)
        
        main_label = ifelse(add_labels, sprintf("bin %s, n=%s", i, sum(f, na.rm=T)), "")
        
        plot(le1[f,1][!(f_weak|f_strong)], le1[f,2][!(f_weak|f_strong)], pch=19, cex=0.5, cex.axis = 2, xlim=c(-20,-14), ylim=c(-20,-14), #cex.main = 0.5, 
             main = main_label, xlab="WT", ylab="DNMT TKO", col = background_color, xaxt = "n", yaxt = "n")
        
        axis(1, labels = FALSE)
        axis(2, labels = FALSE)
        
			if (add_labels) {
				text(x=-20, y=-14.2, labels=paste0(up, "%"), pos=4, cex=2)
				text(x=-14, y=-20, labels=paste0(down, "%"), pos=2, cex=2)
                

			}
        
        points(le1[f,1][f_weak], le1[f,2][f_weak], pch=19, cex=0.5, col = weak_color)
        points(le1[f,1][f_strong], le1[f,2][f_strong], pch=19, cex=0.5, col = strong_color)
        points(le1[f,1][f_up & f_germline], le1[f,2][f_up & f_germline], pch=19, cex=0.5, col = "deepskyblue2")
        
		abline(a=0,b=1)
		abline(a=1,b=1, lty = 2)
		abline(a=-1,b=1, lty = 2)
        f_up = (le1[f,4]-le1[f,1]) > 1
        f_down = (le1[f,4]-le1[f,1]) < -1
		up = round(100*sum(f_up,na.rm=T)/sum(f,na.rm=T), 1)
		down = round(100*sum(f_down,na.rm=T)/sum(f,na.rm=T), 1)
        
        plot(le1[f,1][!(f_weak|f_strong)], le1[f,4][!(f_weak|f_strong)], pch=19, cex=0.5, cex.axis = 2, xlim=c(-20,-14), ylim=c(-20,-14), #cex.main = 0.5, 
             main = main_label, xlab="WT", ylab="DNMT TKO", col = background_color, xaxt = "n", yaxt = "n")
        
        axis(1, labels = FALSE)
        axis(2, labels = FALSE)

			if (add_labels) {
				text(x=-20, y=-14.2, labels=paste0(up, "%"), pos=4, cex=2)
				text(x=-14, y=-20, labels=paste0(down, "%"), pos=2, cex=2)
                
			}
        
        points(le1[f,1][f_weak], le1[f,4][f_weak], pch=19, cex=0.5, col = weak_color)
        points(le1[f,1][f_strong], le1[f,4][f_strong], pch=19, cex=0.5, col = strong_color)
        points(le1[f,1][f_up & f_germline], le1[f,4][f_up & f_germline], pch=19, cex=0.5, col = "deepskyblue2")
        
		abline(a=0,b=1)
		abline(a=1,b=1, lty = 2)
		abline(a=-1,b=1, lty = 2)
	}
	})
	
}
