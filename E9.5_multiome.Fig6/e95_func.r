

# ---- used by prep/E9.5_atac_peaks.ipynb ----

join_norm_peaks.e95 <- function (T_q = 0.97, eps = 3, state = "marginal") 
{
    norm_suf = "marginal_regcap"
    peaks_host = comp_peaks_over_punc_region("embdnmt.brain_e95_host.marginal",
                                             "embdnmt.brain_e95_host.marginal_regcap",
                                             peak_ext = 140, eps = eps, T_q = T_q)
    peaks_dko = comp_peaks_over_punc_region("embdnmt.brain_e95_dko.marginal",
                                            "embdnmt.brain_e95_dko.marginal_regcap",
                                            peak_ext = 140, eps = eps, T_q = T_q)
    peaks = gintervals.union(peaks_host, peaks_dko)
    return(peaks)
}

# As collect_cov_r.e95, but reads npeaks.RDS itself and writes relative to the working directory.
collect_cov_r_npeaks.e95 <- function () 
{
    npeaks = readRDS(in_path("./revision_e95/data/output/npeaks.RDS"))
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

    saveRDS(cov_stat, file = out_path("./revision_e95/data/output/multi_type_cov.RDS"))
    saveRDS(cov_peak_chr, file = out_path("./revision_e95/data/output/cov_peak_chr.RDS"))
    saveRDS(cov_chr, file = out_path("./revision_e95/data/output/cov_chr.RDS"))
    rownames(cov_stat) = cov_stat$nm
}

# ---- used by E9.5_multiome.Fig6/Fig6.ipynb (Fig. 6d) ----

calc_chisq_pval_per_cell_type <- function(umi_df_long, UMI_thresh = 10, num_proc = 40){
    doMC::registerDoMC(num_proc)
    chisq_results <- umi_df_long %>% 
        filter(!(ko + host < UMI_thresh),
               rowSums(. < 0) == 0) %>%
        plyr::adply(1, function(x) {
            chisq <- chisq.test(
                matrix(c(x$ko, x$tot_min_ko, x$host, x$tot_min_host), byrow = TRUE, nrow = 2)
            )
            tibble(pval = chisq$p.value, statistic = chisq$statistic)
        }, .parallel = TRUE) %>%
        group_by(cell_type) %>% 
        mutate(qval = p.adjust(pval, method = "fdr")) %>% 
        ungroup()
    
    return(chisq_results)
    }

plot_brain_expr <- function(gene, cell_types, plot_sig = TRUE, lower_y = log2(2e-5), upper_y = -7.5){

  highlight_gn <- as.character(gene)
  cell_types_chr <- trimws(as.character(cell_types))
    cell_types_chr = cell_types_chr[cell_types_chr %in% colnames(mko_l)]

  plot_df <- data.frame(
    ko   = mko_l[highlight_gn, cell_types_chr],
    host = mhost_l[highlight_gn, cell_types_chr]
  ) %>%
    rownames_to_column("cell_type") %>%
    mutate(cell_type_chr = trimws(as.character(cell_type)),
           x_id = match(cell_type_chr, cell_types_chr)) %>%
    filter(!is.na(x_id))

  df_long <- plot_df %>%
    pivot_longer(c(ko, host), names_to = "group", values_to = "value") %>%
    mutate(
      group = factor(group, levels = c("host", "ko"))
    )

  # geometry parameters
  dodge_w <- 0.35
  bar_w   <- 0.3
  min_h   <- 0.08

  df_rect <- df_long %>%
    mutate(
      x_center = x_id + if_else(group == "host", -dodge_w/2, dodge_w/2),
      xmin = x_center - bar_w/2,
      xmax = x_center + bar_w/2,
      ymin = lower_y,
      ymax = if_else(value < lower_y, lower_y + min_h, value)
    )

    if (plot_sig){
          df_y <- df_rect %>%
            group_by(x_id) %>%
            summarise(y_pair = max(ymax, na.rm = TRUE), .groups = "drop")

          df_p <- chisq %>%
            mutate(
              gene = trimws(as.character(gene)),
              cell_type_chr = trimws(as.character(cell_type)),
              x_id = match(cell_type_chr, cell_types_chr),
              p_lab = case_when(
                is.na(qval)   ~ "NA",
                qval < 0.001  ~ "***",
                qval < 0.01   ~ "**",
                qval < 0.05   ~ "*",
                TRUE          ~ "NS"
              )
            ) %>%
            filter(gene == highlight_gn, !is.na(x_id))
        
          # ---- bracket coords ----
  tick_h  <- 0.12
  text_up <- 0.4

  df_br <- df_p %>%
    left_join(df_y, by = "x_id") %>%
    mutate(
      y_pair = if_else(is.na(y_pair), lower_y + min_h, y_pair),
      gap_up = 0.22,
      y = y_pair + gap_up,
      y = pmin(pmax(y, lower_y + 0.25), upper_y - 0.35),
      x1 = x_id - dodge_w/2 - bar_w/2,
      x2 = x_id + dodge_w/2 + bar_w/2,
      xm = x_id,
      y_text = pmin(y + text_up, upper_y - 0.10)
    )
        
    }

    message(sprintf("%s: %s genes too high.", gene, sum(df_rect$ymax > upper_y)))

  p <- ggplot(df_rect) +
    geom_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = cell_type_chr)) +
    geom_hline(yintercept = log2(2e-5)) + 
    scale_fill_manual(values = mod_e95$colors_e95, drop = FALSE) +
    scale_x_continuous(
      limits = c(0.5, length(cell_types_chr) + 0.5),
      breaks = seq_along(cell_types_chr),
      labels = cell_types_chr,
      expand = c(0, 0)
    ) +
    scale_y_continuous(
      limits = c(lower_y, upper_y),
      expand = c(0, 0)
    ) +
    labs(x = NULL, y = NULL, title = highlight_gn) +
    theme_bw(base_size = 10) +
    theme(
      plot.title = element_text(size = 12),
      axis.text.x = element_blank(),
      axis.ticks  = element_blank(),
      panel.grid.major.x = element_blank(),
      legend.position = "none"
    )
    
    if (plot_sig){
        p = p +
            geom_segment(
              data = df_br, inherit.aes = FALSE,
              aes(x = x1, xend = x2, y = y, yend = y)
            ) +
            geom_segment(
              data = df_br, inherit.aes = FALSE,
              aes(x = x1, xend = x1, y = y, yend = y - tick_h)
            ) +
            geom_segment(
              data = df_br, inherit.aes = FALSE,
              aes(x = x2, xend = x2, y = y, yend = y - tick_h)
            ) +
            geom_text(
              data = df_br, inherit.aes = FALSE,
              aes(x = xm, y = y_text, label = p_lab),
              size = 3
            )
    }
    
    return(p)
}

# ---- used by E9.5_multiome.Fig6/Fig6.ipynb (Fig. 6b, c) ----

darken_hex <- function(hex, factor = 0.82) {
  rgb_mat <- grDevices::col2rgb(hex) / 255
  rgb_new <- pmax(pmin(rgb_mat * factor, 1), 0)
  grDevices::rgb(rgb_new[1,], rgb_new[2,], rgb_new[3,])
}

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
