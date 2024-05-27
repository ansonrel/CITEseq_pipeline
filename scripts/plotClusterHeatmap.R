#' @rdname plotClusterHeatmap
#' @title Plot cluster heatmap
#' 
#' @description 
#' Plots heatmaps summarizing a clustering and/or metaclustering of interest.
#'
#' @param x a \code{\link[SingleCellExperiment]{SingleCellExperiment}}.
#' @param hm2 character string. Specifies the right-hand side heatmap. 
#'   One of: \itemize{
#'   \item{\code{"abundances"}: cluster frequencies across samples}
#'   \item{\code{"state_markers"}: median cell state marker expressions 
#'     across clusters (analogous to the left-hand side heatmap)}
#'   \item{a character string/vector corresponding to one/multiple marker(s): 
#'     median marker expressions across samples and clusters}}
#' @param k 
#'   character string. Specifies the clustering 
#'   across which median marker expressions should be computed.
#' @param m 
#'   character string. Specifies the metaclustering to be shown. 
#'   (This is for display only and will not effect any computations!) 
#' @param fun
#'   character string specifying the function to use as summary statistic.
#' @param cluster_anno 
#'   logical. Specifies if clusters should be annotated.
#' @param split_by 
#'   character string. Must corresponds to a column name of \code{rowData(x)}. 
#'   If specified, the data will be subset according to this variable, 
#'   and multiple heatmaps will be drawn.
#' @param scale 
#'   logical. Specifies whether scaled values should be plotted.
#'   (see below for details)
#' @param draw_dend 
#'   logical. Specifies if the row dendrogram should be drawn.
#' @param draw_freqs 
#'   logical. Specifyies whether to display cell counts and proportions.
#' @param palette 
#'   character vector of colors to interpolate. 
#' 
#' @details Scaled values corresponds to cofactor arcsinh-transformed 
#' expression values scaled between 0 and 1 using 1% and 99% percentiles as 
#' boundaries. Hierarchical clustering is performed on the unscaled data.
#' 
#' In its 1st panel, \code{plotClusterHeatmap} will display
#' median (scaled, arcsinh-transformed) cell-type marker expressions (across all samples).
#' Depending on argument \code{hm2}, the 2nd panel will contain one of:
#' \itemize{
#' \item{relataive cluster abundances by sample}
#' \item{median (scaled, arcsinh-transformed) 
#'   cell-state marker expressions (across all samples)}
#' \item{median (scaled, arcsinh-transformed) 
#'   cell-state marker expressions by sample}
#' }
#' 
#' @return a \code{\link{HeatmapList-class}} object.
#' 
#' @author Helena Lucia Crowell \email{helena.crowell@@uzh.ch}
#' 
#' @references 
#' Nowicka M, Krieg C, Weber LM et al. 
#' CyTOF workflow: Differential discovery in 
#' high-throughput high-dimensional cytometry datasets.
#' \emph{F1000Research} 2017, 6:748 (doi: 10.12688/f1000research.11622.1)
#' 
#' @examples
#' # construct SCE & run clustering
#' data(PBMC_fs, PBMC_panel, PBMC_md)
#' sce <- prepData(PBMC_fs, PBMC_panel, PBMC_md)
#' sce <- cluster(sce)
#' 
#' plotClusterHeatmap(sce, hm2="abundances")
#' plotClusterHeatmap(sce, hm2="abundances", draw_freqs=TRUE)
#' plotClusterHeatmap(sce, hm2="state_markers", k="meta16", split_by='condition')
#' plotClusterHeatmap(sce, hm2="pS6", k="meta12", m="meta8")
#' plotClusterHeatmap(sce, hm2="abundances", scale=FALSE, draw_freqs=TRUE)
#' 
#' @import ComplexHeatmap
#' @importFrom grDevices colorRampPalette
#' @importFrom RColorBrewer brewer.pal
#' @importFrom S4Vectors metadata
#' @importFrom stats dist hclust
#' @importFrom SummarizedExperiment assay
#' @export

plotClusterHeatmap <- function(x, 
                               hm2=NULL, 
                               #k="meta20", 
                               m=NULL, fun=c("median", "mean"), 
                               cluster_anno=TRUE, split_by=NULL, scale=TRUE, 
                               draw_dend=TRUE, draw_freqs=FALSE, 
                               palette=rev(brewer.pal(11, "RdYlBu"))) {
  
  # check validity of input arguments
  # .check_sce(x)
  # fun <- match.arg(fun)
  # k <- .check_validity_of_k(x, k)
  # .check_cd_factor(x, split_by)
  # u <- c("abundances", "state_markers", rownames(x))
  # if (!is.null(hm2)) stopifnot(hm2 %in% u)
  
  # hierarchical clustering on cell-type marker medians by cluster_id
  #nk <- nlevels(x$cluster <- cluster_ids(x, k))
  nk <- nlevels(x$cluster_id)
  
  
  .split_cells <- function(x, by) {
    stopifnot(is.character(by), by %in% colnames(colData(x)))
    cd <- data.frame(colData(x))
    dt <- data.table(cd, i = seq_len(ncol(x)))
    dt_split <- split(dt, by = by, sorted = TRUE, flatten = FALSE)
    map_depth(dt_split, length(by), "i")
  }
  
  .agg <- function(x, by, fun = c("median", "mean", "sum")) {
    fun <- switch(match.arg(fun), 
                  median = rowMedians, mean = rowMeans, sum = rowSums)
    cs <- .split_cells(x, by)
    pb <- map_depth(cs, -1, function(i) {
      if (length(i) == 0) return(numeric(nrow(x)))
      fun(assay(x, "exprs")[, i, drop = FALSE])
    })
    map_depth(pb, -2, function(u) as.matrix(data.frame(
      u, row.names = rownames(x), check.names = FALSE)))
  }
  
  
  ms_by_k <- t(.agg(x, "cluster_id", fun))
  
  #d <- dist(ms_by_k[, type_markers(x)])
  d <- dist(ms_by_k)
  row_clustering <- hclust(d, method="average")
  
  # clustering row annotation
  
  .cluster_cols <- c(
    "#DC050C", "#FB8072", "#1965B0", "#7BAFDE", "#882E72",
    "#B17BA6", "#FF7F00", "#FDB462", "#E7298A", "#E78AC3",
    "#33A02C", "#B2DF8A", "#55A1B1", "#8DD3C7", "#A6761D",
    "#E6AB02", "#7570B3", "#BEAED4", "#666666", "#999999",
    "#aa8282", "#d4b7b7", "#8600bf", "#ba5ce3", "#808000",
    "#aeae5c", "#1e90ff", "#00bfff", "#56ff0d", "#ffff00")
  
  
  if (cluster_anno) {
    anno <- levels(x$cluster_id)
    if (nk > 30) {
      
      cols <- colorRampPalette(.cluster_cols)(nk)
    } else {
      cols <- .cluster_cols[seq_len(nk)]
    }
    cols <- setNames(cols, anno)
    
    .row_anno <- function(anno, cols, name, clustering, dend) {
      Heatmap(matrix=anno, col=cols, name=name, 
              rect_gp=gpar(col="white"), width=unit(.4, "cm"),
              cluster_rows=clustering, cluster_columns=FALSE,
              show_row_dend=dend, row_dend_reorder=FALSE)
    }
    
    cluster_anno <- .row_anno(anno, cols, 
                              "cluster_id", row_clustering, draw_dend)
  }
  # merging row annotation
  if (!is.null(m)) {
    #.check_validity_of_k(x, m)
    idx <- match(seq_len(nk), cluster_codes(x)[, k])
    anno <- factor(cluster_codes(x)[, m][idx])
    if (nlevels(anno) > 30) {
      cols <- colorRampPalette(.cluster_cols)(nlevels(anno))
    } else {
      cols <- .cluster_cols[seq_len(nlevels(anno))]
    }
    cols <- setNames(cols, levels(anno))
    merging_anno <- .row_anno(anno, cols,
                              "merging_id", row_clustering, draw_dend)
  }

  # split cell indices by colData factor
  many <- !is.null(split_by)
  cs <- seq_len(ncol(x))
  if (many) groups <- split(cs, x[[split_by]]) else groups <- list(cs)   
  
  # optionally scale expression matrix
  
  .scale_exprs <- function(x, margin = 1) {
    if (!is(x, "matrix")) x <- as.matrix(x)
    qs <- c(rowQuantiles, colQuantiles)[[margin]]
    qs <- qs(x, probs = c(.01, .99))
    x <- switch(margin,
                "1" = (x - qs[, 1]) / (qs[, 2] - qs[, 1]),
                "2" = t((t(x) - qs[, 1]) / (qs[, 2] - qs[, 1])))
    x[x < 0 | is.na(x)] <- 0
    x[x > 1] <- 1
    return(x)
  }
  
  if (scale) assay(x, "exprs")  <- 
    .scale_exprs(assay(x, "exprs") , 1)
  hm_cols <- colorRampPalette(palette)(100)
  
  hms <- lapply(seq_along(groups), function(i) {
    idx <- groups[[i]]
    cs_by_k <- split(idx, x$cluster_id[idx])
    # left-hand side heatmap:
    # median cell-type marker expressions across clusters
    if (!many) {
      if (scale) {
        hm1_es <- t(.agg(x, "cluster_id", fun))
      } else {
        hm1_es <- ms_by_k
      }
    } else {
      hm1_es <- t(.agg(x[, idx], "cluster_id", fun))
    }
    hm1 <- Heatmap(
      #matrix=hm1_es[, type_markers(x)],
      matrix=hm1_es, 
      col=hm_cols, name="expression", 
      column_names_gp=gpar(fontsize=8),
      rect_gp=gpar(col='white'), na_col="lightgrey", 
      cluster_rows=row_clustering, cluster_columns=FALSE,
      show_row_dend=draw_dend, column_title=names(groups)[i][many])
    
    # cluster frequencies
    freq_bars <- freq_anno <- NULL
    if (draw_freqs) {
      fq <- round(tabulate(x$cluster_id[idx]) / length(idx) * 100, 2)
      freq_bars <- rowAnnotation(
        "Frequency [%]"=row_anno_barplot(
          x=fq, axis=TRUE, border=FALSE, bar_with=.8, 
          gp=gpar(fill="grey50", col="white")), width=unit(2, "cm"))
      labs <- paste0(levels(x$cluster_id), " (", fq, "%)")
      freq_anno <- rowAnnotation(
        text=row_anno_text(labs), 
        width=max_text_width(labs))
    }
    
    # combine row annotations, heatmap, 
    # and frequency bars & labels
    p <- hm1 + freq_bars + freq_anno
    if (is(cluster_anno, "Heatmap")) 
      p <- cluster_anno + p
    if (exists("merging_anno")) 
      p <- merging_anno + p
    
    # right-hand side heatmap
    if (!is.null(hm2)) {
      if (hm2 == "abundances") {
        # cluster frequencies across samples
        cs <- table(x$cluster_id[idx], x$sample_id[idx])
        fq <- as.matrix(unclass(prop.table(cs, 2)))
        fq <- fq[, !is.na(colSums(fq)), drop = FALSE]
        p <- p + Heatmap(matrix=fq, name="frequency", 
                         na_col="lightgrey", rect_gp=gpar(col="white"), 
                         show_row_names=FALSE, column_names_gp=gpar(fontsize=8), 
                         cluster_rows=row_clustering, cluster_columns=FALSE)
      } else if (hm2 == "state_markers") {
        # median cell state marker expressions across clusters
        p <- p + Heatmap(col=hm_cols, na_col="lightgrey", 
                         matrix=hm1_es[, state_markers(x)], 
                         rect_gp=gpar(col='white'), show_heatmap_legend=FALSE, 
                         cluster_rows=row_clustering, cluster_columns=FALSE,
                         column_names_gp=gpar(fontsize=8))
      } else {
        for (ch in hm2) {
          # median marker expression across samples & clusters
          ms <- .agg(x[ch, idx], c("cluster_id", "sample_id"), fun)
          ms <- do.call("rbind", ms)
          rownames(ms) <- levels(x$cluster_id)
          p <- p + Heatmap(matrix=ms, col=hm_cols, 
                           na_col="lightgrey", rect_gp=gpar(col='white'),
                           show_heatmap_legend=FALSE, show_row_names=FALSE,
                           cluster_rows=row_clustering, cluster_columns=FALSE,
                           column_title=ch, column_names_gp=gpar(fontsize=8))
        }
      }
    }
    return(p)
  })
  for (i in seq_along(hms)) 
    draw(hms[[i]])
  invisible(hms)
}