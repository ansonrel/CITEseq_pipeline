
# custom plotDiffHeatmap that groups results by cell-types and allows to sort by p_adj.loc
# adapted from https://github.com/HelenaLC/CATALYST/blob/master/R/plotDiffHeatmap.R
# to check the moditications, "CTRL + F ###"
# still needs to have CATALYST installed because the script loads hidden functions (`CATALYST:::.fun`)
# 
# new requirement: 
# - `cluster_id` in colData(x) and the corresponding cluster_codes in metadata(x)
# 

custom_plotDiffHeatmap <- function(x, y, k = NULL,
    top_n = 20, fdr = 0.05, lfc = 1, all = FALSE,
    sort_by = c("padj", "p_adj.loc", "lfc", "none"), ### <--- now allows to sort by p_adj.loc
    y_cols = list(padj = "p_adj", lfc = "logFC", target = "marker_id"),
    assay = "exprs", fun = c("median", "mean", "sum"), 
    normalize = TRUE, col_anno = TRUE, row_anno = TRUE,
    hm_pal = NULL, 
    fdr_pal = c("lightgrey", "lightgreen"),
    lfc_pal = c("blue3", "white", "red3")) {
    
    # check validity of input arguments
    fun <- match.arg(fun)
    sort_by <- match.arg(sort_by)
    args <- as.list(environment())
    
    defs <- as.list(formals("plotDiffHeatmap")$y_cols[-1])
    miss <- !names(defs) %in% names(args$y_cols)
    if (any(miss)) y_cols <- args$y_cols <- 
        c(args$y_cols, defs[miss])[names(defs)]

    ### ---------->
    # removes checks. Causes problems because of new argument (sort_by = "p_adj.loc")
    #.check_args_plotDiffHeatmap(args)
    #stopifnot(y_cols[[sort_by]] %in% names(y))
    #y_cols <- y_cols[y_cols %in% names(y)]
    ### <-----------
        
    # guess clustering to use
    if (is.null(k)) {
        kids <- levels(y$cluster_id)
        same <- vapply(cluster_codes(x), function(u) 
            identical(levels(u), kids), logical(1))
        if (!any(same)) 
            stop("Couldn't match any clustering",
                " in input data 'x' with results in 'y'.")
        k <- names(cluster_codes(x))[same][1]
    } else {
        k <- CATALYST:::.check_k(x, k)
    }
    x$cluster_id <- cluster_ids(x, k)
    
    # get feature column
    y <- data.frame(y, check.names = FALSE)
    y <- dplyr::mutate_if(y, is.factor, as.character)
    if (any(rownames(x) %in% unlist(y))) {
        features <- intersect(rownames(x), y[[y_cols$target]])
        if (length(features) == 0)
            stop("Couldn't match features between",
                " results 'y' and input data 'x'.")
        i <- y[[y_cols$target]] %in% rownames(x)
        type <- "ds"
    } else {
        i <- TRUE
        type <- "da"
    }

    # rename relevant result variables
    y <- dplyr::rename(y, 
        target = y_cols$target,
        padj = y_cols$padj, 
        lfc = y_cols$lfc)
    
    # filter results
    i <- i & !is.na(y$padj) & y$cluster_id %in% levels(x$cluster_id)
    if (!all) {
        i <- i & y$padj < fdr
        if (!is.null(y$lfc))
            i <- i & abs(y$lfc) > lfc
    }
    y <- y[i, , drop = FALSE]

    if (nrow(y) == 0)
        stop("No results remaining;",
            " perhaps 'x' or 'y' has been filtered,",
            " or features couldn't be matched.")
    
    # get clusters/cluster-marker combinations to plot
    if (sort_by != "none") {
        o <- order(abs(y[[sort_by]]), 
            decreasing = (sort_by == "lfc"))
        y <- y[o, , drop = FALSE]
    }
    if (top_n > nrow(y)) 
        top_n <- nrow(y)
    top <- y[seq_len(top_n), ]
    
    # column annotation of non-numeric cell metadata variables
    if (!isFALSE(col_anno)) {
        top_anno <- CATALYST:::.anno_factors(x, levels(x$sample_id), col_anno, "column")
    } else top_anno <- NULL
    
    if (is.null(hm_pal)) hm_pal <- rev(RColorBrewer::brewer.pal(11, 
        ifelse(type == "ds", "RdYlBu", "RdBu")))
    
    # row annotation: significant = (adj. p-values <= th)
    if (row_anno) {
        s <- factor(
            ifelse(top$padj < fdr, "yes", "no"), 
            levels = c("no", "yes"))
        if (!is.null(top$lfc)) {
            lfc_lims <- range(top$lfc, na.rm = TRUE)
            if (all(lfc_lims > 0)) {
                lfc_brks <- c(0, lfc_lims[2])
                lfc_pal <- lfc_pal[-1]
            } else if (all(lfc_lims < 0)) {
                lfc_brks <- c(lfc_lims[1], 0)
                lfc_pal <- lfc_pal[-3]
            } else lfc_brks <- c(lfc_lims[1], 0, lfc_lims[2])
            lfc_anno <- top$lfc
            anno_cols <- list(logFC = circlize::colorRamp2(lfc_brks, lfc_pal))
        } else {
            lfc_anno <- NULL
            anno_cols <- list()
        }


        names(fdr_pal) <- levels(s)
        anno_cols$significant <- fdr_pal

        right_anno <- ComplexHeatmap::rowAnnotation(
            logFC = lfc_anno,
            significant = s,
            "foo" = ComplexHeatmap::row_anno_text(
                scales::scientific(top$padj, 2),
                gp = grid::gpar(fontsize = 8)),
            col = anno_cols,
            gp = grid::gpar(col = "white"),
            show_annotation_name = FALSE,
            simple_anno_size = grid::unit(4, "mm"))
        left_color = structure(1:length(unique(top$cluster_id)), 
                    names = unique(top$cluster_id)) # black, red, green, blue
        ### ------------->
        # Add left annotation heatmap for cell-types. Needs `cluster_id in colData(x)
        left_anno <- ComplexHeatmap::rowAnnotation(
            cluster_id = top$cluster_id,
            col = list(cluster_id = left_color),
            show_annotation_name = FALSE,
            simple_anno_size = grid::unit(4, "mm"))
        ### <-----------
    } else right_anno <- NULL

    switch(type, 
        # relative cluster abundances by sample
        da = {
            ns <- table(x$cluster_id, x$sample_id)
            fq <- prop.table(ns, 2)
            fq <- fq[top$cluster_id, ]
            y <- as.matrix(unclass(fq))
            if (normalize) y <- CATALYST:::.z_normalize(asin(sqrt(y)))
            ComplexHeatmap::Heatmap(
                matrix = y, 
                name = paste0("normalized\n"[normalize], "frequency"),
                col = hm_pal,
                na_col = "lightgrey", 
                rect_gp = grid::gpar(col = "white"),  
                cluster_rows = FALSE,
                cluster_columns = FALSE,
                row_names_side = "left",
                top_annotation = top_anno,
                right_annotation = right_anno)
        },
        # median state-marker expression by sample
        ds = {
            y <- assay(x, assay)
            cs <- CATALYST:::.split_cells(x, c("cluster_id", "sample_id"))
            z <- t(mapply(function(k, g)
                vapply(cs[[k]], function(cs) {
                    if (length(cs) == 0) return(NA)
                    get(fun)(y[g, cs, drop = FALSE])
                }, numeric(1)),
                k = top$cluster_id, 
                g = top$target))
            rownames(z) <-  top$target
            if (normalize) z <- CATALYST:::.z_normalize(z) 
            ComplexHeatmap::Heatmap(
                matrix = z,
                name = paste0("z-normalized\n"[normalize], "expression"),
                col = hm_pal,
                cluster_rows = FALSE,
                row_split = top$cluster_id,
                cluster_columns = FALSE,
                top_annotation = top_anno,
                row_names_side = "left",
                rect_gp = grid::gpar(col = "white"),
                right_annotation = right_anno,
                left_annotation = left_anno, ### <----- left annotation
                heatmap_legend_param = list(title_gp = grid::gpar(
                    fontsize = 10, fontface = "bold", lineheight = 0.8)))
    }) 
}
