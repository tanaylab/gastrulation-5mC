# Functions used by ExtData_RNA.ipynb.

# ---- used by ExtData_RNA.ipynb ----

plot_marker_htmp <- function(name, gns, width = 7.5, height = 7.5){
    
    annotation_col = RNA_mod[[name]]$mc_to_type %>% enframe(value = "cell_type") %>% column_to_rownames("name")
    cts_to_keep = setdiff(names(annotation_colors$cell_type), c('Visceral endoderm','ExE visceral endoderm','Visceral endoderm_anterior'))
    annotation_col = annotation_col %>% filter(cell_type %in% cts_to_keep)
    annotation_col = annotation_col[order(match(annotation_col$cell_type, names(annotation_colors$cell_type))),,F]
    
    gns_to_plot = gns[gns %in% rownames(RNA_mod[[name]]$legc)]
    mat = RNA_mod[[name]]$legc[gns_to_plot, rownames(annotation_col)]
    fc_mat = t(apply(mat, 1, function(x) {x - median(x)}))

    # Set zero to always be white
    range_data <- range(fc_mat)
    min_data <- -2 # range_data[1]
    max_data <- 2 # range_data[2]
    breaks <- c(seq(min_data, 0, length.out = 51), seq(0, max_data, length.out = 50)[-1])

    custom_palette <- colorRampPalette(c("black", "grey30", "grey60", "grey90", "white", "deepskyblue1", "deepskyblue2", "deepskyblue3", "deepskyblue4"))(length(breaks) - 1)
    
    show_fig(sprintf("./figs/Supp/marker_heatmaps/%s.pdf", name), draw = pheatmap::pheatmap(fc_mat, annotation_col = annotation_col, annotation_colors = annotation_colors, 
                       annotation_legend = F, show_colnames = F, 
                       color = custom_palette, breaks = breaks,
                       w=width, h=height,
                       treeheight_col = 0, cluster_rows = F, cluster_cols = F))
    
}

plot_abs_htmp <- function(name, gns, width = 7.5, height = 7.5){
    
    annotation_col = RNA_mod[[name]]$mc_to_type %>% enframe(value = "cell_type") %>% column_to_rownames("name")
    cts_to_keep = setdiff(names(annotation_colors$cell_type), c('Visceral endoderm','ExE visceral endoderm','Visceral endoderm_anterior'))
    annotation_col = annotation_col %>% filter(cell_type %in% cts_to_keep)
    annotation_col = annotation_col[order(match(annotation_col$cell_type, names(annotation_colors$cell_type))),,F]
    
    gns_to_plot = gns[gns %in% rownames(RNA_mod[[name]]$legc)]
    mat = RNA_mod[[name]]$legc[gns_to_plot, rownames(annotation_col)]

    breaks <- seq(-16, -6, length.out = 1000)

    custom_palette <- colorRampPalette(c("white", "deepskyblue1", "deepskyblue2", "deepskyblue3", "deepskyblue4", "black"))(length(breaks) - 1)
    
    
    p <- show_fig(sprintf("./figs/Supp/imprinted_heatmaps/%s.pdf", name), draw = pheatmap::pheatmap(mat, annotation_col = annotation_col, annotation_colors = annotation_colors, 
                       annotation_legend = F, show_colnames = F, 
                       color = custom_palette, breaks = breaks,
                       w=width, h=height,
                       treeheight_col = 0, cluster_rows = F, cluster_cols = F))
    
    return(invisible(p))
    
}

make_df <- function(expr){

    df = rbind(data.frame("Class" = "Ribosomal", "max" = apply(expr[ribo,], 1, max)),
               data.frame("Class" = "Key TFs", "max" = apply(expr[markers,], 1, max)),
               data.frame("Class" = "Germline genes,\nDahlet et al. 2020", "max" = apply(expr[weber_f,], 1, max))) %>%
         rownames_to_column("gene")

    return(df)
}
