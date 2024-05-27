

getUMAPclusters <- function(d, coldata = "cluster_id", meta = NULL, dr = "UMAP"){
  
coerce_to_factor <- function(x, level.limit) {
  if (!is.null(x)) {
    x <- as.factor(x)
    if (nlevels(x) > level.limit) {
      stop(sprintf("more than %i levels for '%s'", level.limit))
    }
  }
  x
}

red_dim <- reducedDim(d, dr)
colnames(red_dim) <- NULL 
df_to_plot <- data.frame(red_dim)

if(!is.null(meta)) {
  d$k <- factor(cluster_ids(d, meta))
  text_out <- retrieveCellInfo(d, "k", search="colData")
} else {
  text_out <- retrieveCellInfo(d, coldata, search="colData")
}

text_out$val <- coerce_to_factor(text_out$val, level.limit=Inf)
#cluster name locations
by_text_x <- vapply(split(df_to_plot$X1, text_out$val), median, na.rm = TRUE, FUN.VALUE=0)
by_text_y <- vapply(split(df_to_plot$X2, text_out$val), median, na.rm = TRUE, FUN.VALUE=0)

return(list(x = by_text_x, y = by_text_y))
}