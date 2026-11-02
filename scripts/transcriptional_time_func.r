# Functions used by prep/transcriptional_time.ipynb.

# ---- used by prep/transcriptional_time.ipynb ----

get_query_and_atlas_time_dist = function(query_cls_md, atlas_time, graph_id) {

  # query_cls_md is a named vector with query_cls as names and the embryo annotation or sth similar as value
  # atlas_time is a named vector giving the atlas_time per for each atlas cell

  cgraph = scdb_cgraph(graph_id)

  query_cls = intersect(names(query_cls_md), cgraph@nodes)
  atlas_cls = intersect(names(atlas_time), cgraph@nodes)
  atlas_time = atlas_time[atlas_cls]
  query_cls_md = query_cls_md[query_cls]
    
  cell_names = c(1:length(cgraph@nodes))
  names(cell_names) = cgraph@nodes

  graph = cgraph@edges
  graph$mc1 = as.factor(graph$mc1)
  graph$mc2 = as.factor(graph$mc2)
  levels(graph$mc1) = cell_names[levels(graph$mc1)]
  levels(graph$mc2) = cell_names[levels(graph$mc2)]

  knn_mat = sparseMatrix(as.numeric(graph$mc1), as.numeric(graph$mc2),x = graph$w)
  colnames(knn_mat) = cgraph@nodes
  rownames(knn_mat) = cgraph@nodes

  knn_mat_f = knn_mat[atlas_cls,atlas_cls]
  for(i in unique(atlas_time)) {
    f = atlas_time == i
    knn_mat_f[f,f] = 0
  }

  time_match = apply(knn_mat_f,1,which.max)
  time_match = atlas_time[time_match]

  q_time_match_tmp = table(atlas_time, time_match)
  q_time_match = matrix(0,nrow = length(unique(atlas_time)),ncol = length(unique(atlas_time)))
  rownames(q_time_match) = sort(unique(atlas_time))
  colnames(q_time_match) = sort(unique(atlas_time))
  q_time_match[rownames(q_time_match_tmp),colnames(q_time_match_tmp)] = q_time_match_tmp

  query_time_match = apply(knn_mat[query_cls, atlas_cls], 1, which.max)
  query_time_dist = table(query_cls_md, atlas_time[query_time_match])
  tmp = matrix(0,nrow = nrow(query_time_dist),ncol = length(unique(atlas_time)))
  rownames(tmp) = rownames(query_time_dist)
  colnames(tmp) = sort(unique(atlas_time))
  tmp[rownames(query_time_dist),colnames(query_time_dist)] = query_time_dist
  query_time_dist = tmp

  return(list(atlas_time_dist = q_time_match, query_time_dist = query_time_dist))
}

tgs_cor_graph = function(x, knn, k_expand, k_alpha, k_beta)
{
	x_knn = tgs_cor_knn(x, y=x, knn=knn*k_expand)
	gr = tgs_graph(x_knn, knn, k_expand, k_beta)
	return(gr)
}

tgs_cor_knn <- function(x, y, knn, pairwise.complete.obs = FALSE, spearman = FALSE, threshold = 0) {
    if (missing(x) || missing(knn)) {
        stop("Usage: tgs_cor_knn(x, y, knn, pairwise.complete.obs = FALSE, spearman = FALSE, threshold = 0)", call. = FALSE)
    }

    if (is.null(y)) {
        .Call("tgs_cor_knn", x, knn, pairwise.complete.obs, spearman, threshold, new.env(parent = parent.frame()))
    } else {
        .Call("tgs_cross_cor_knn", x, y, knn, pairwise.complete.obs, spearman, threshold, new.env(parent = parent.frame()))
    }
}
