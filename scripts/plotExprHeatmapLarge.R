#' @importFrom ComplexHeatmap Heatmap
#' @importFrom circlize colorRamp2
#' @importFrom dplyr select select_if
#' @importFrom grid gpar grid.text
#' @importFrom grDevices colorRampPalette
#' @importFrom magrittr set_rownames
#' @importFrom RColorBrewer brewer.pal
#' @importFrom S4Vectors metadata
#' @importFrom stats dist hclust
#' @importFrom SummarizedExperiment assay
plotExprHeatmapFull <- function(x, features = NULL, 
                                drow = NULL, 
                            assay = "exprs",
                            cells = 40,
                            k = "cluster_id",
                            annot_cols = c("sample_id", "cluster_id"),
                            split = NULL, 
                            col_anno = TRUE,
                            row_clust = TRUE, col_clust = TRUE, 
                            row_dend = TRUE, col_dend = TRUE, 
                            hm_pal = rev(brewer.pal(11, "RdYlBu")), 
                            distance = c(
                              "euclidean", "maximum", "manhattan", 
                              "canberra", "binary", "minkowski"), 
                            linkage = c(
                              "average", "ward.D", "single", "complete", 
                              "mcquitty", "median", "centroid", "ward.D2")) {
  
  cluster_cols <- c(
    "#DC050C", "#FB8072", "#1965B0", "#7BAFDE", "#882E72",
    "#B17BA6", "#FF7F00", "#FDB462", "#E7298A", "#E78AC3",
    "#33A02C", "#B2DF8A", "#55A1B1", "#8DD3C7", "#A6761D",
    "#E6AB02", "#7570B3", "#BEAED4", "#666666", "#999999",
    "#aa8282", "#d4b7b7", "#8600bf", "#ba5ce3", "#808000",
    "#aeae5c", "#1e90ff", "#00bfff", "#56ff0d", "#ffff00")
  
  # check validity of input arguments
  args <- as.list(environment())
  #.check_args_plotExprHeatmap(args)
  distance <- match.arg(distance)
  linkage <- match.arg(linkage)
  
  # subset features of interest
  x <- x[features,]
  
  # use specified cluster id
  x$cluster_id <- x[[k]] 
  
  # Choose cells per cluster
  
  if (is.null(cells)) {
  # use all cells
  cs <- TRUE 
} else {
  # if (is.null(x$sample_id))
  #   stop("colData column sample_id not found,\n ", 
  #        " but is required to downsample cells.")
  stopifnot(
    is.numeric(cells), length(cells) == 1,
    as.integer(cells) == cells, cells > 0)
  
  # split cell indices by cluster
  cs <- split(seq_len(ncol(x)), x$cluster_id)
  # sample at most 'n' cells per sample
  cs <- unlist(lapply(cs, function(u)
    sample(u, min(cells, length(u)))))
}

  #column anno
  col_fun <- 
  if (!isFALSE(col_anno)) {
    top_anno <- anno_clusters(x, k_pal = cluster_cols, idx = cs) 
  } else top_anno <- NULL
  
  
  #split by columns 
  if(!is.null(split)) {
    split <- split[cs]
  }
  z <- as.matrix(assay(x[,cs], assay))
  
  #row anno
  side_anno <- anno_rows(x, k_pal = cluster_cols, genesdf = drow)

  #draw heatmap  
  Heatmap(
    matrix = z,
    name = assay,
    col = colorRamp2(
      seq(min(z), max(z), l = n <- 100),
      colorRampPalette(hm_pal)(n)),
    cluster_rows = row_clust,
    cluster_columns = col_clust,
    show_row_dend = row_dend,
    show_column_dend = col_dend,
    clustering_distance_rows = distance,
    clustering_method_rows = linkage,
    clustering_distance_columns = distance,
    clustering_method_columns = linkage,
    show_row_names = TRUE,
    show_column_names = FALSE,
    top_annotation = top_anno,
    column_title = "cells",
    column_title_side = "bottom",
    column_split = split,
    left_annotation = side_anno
    #right_annotation = right_anno,
    #heatmap_legend_param = lgd_aes
  )

}


anno_clusters <- function(x, k_pal, idx) {
  kids <- levels(x$cluster_id)
  nk <- length(kids)
  if (nk > length(k_pal)) k_pal <- colorRampPalette(k_pal)(nk)
  k_pal <- k_pal[seq_len(nk)]
  names(k_pal) <- kids
  #df <- data.frame(cluster_id = kids)
  col <- list(cluster_id = k_pal)
  
  df <- data.frame(cluster_id = colData(x)[idx,"cluster_id"])
  HeatmapAnnotation(which = "column", 
                    df = df,
                    col = col)
  
}

anno_rows <- function(x, k_pal, genesdf) {
  kids <- levels(x$cluster_id)
  nk <- length(kids)
  if (nk > length(k_pal)) k_pal <- colorRampPalette(k_pal)(nk)
  k_pal <- k_pal[seq_len(nk)]
  names(k_pal) <- kids
  col <- list(cluster_id = k_pal)
  
  row_ha <- HeatmapAnnotation(which = "row", df = genesdf, col = col)
}
