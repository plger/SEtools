#' crossHm
#'
#' Plot a multi-panel heatmap from a list of
#' \code{\link[SummarizedExperiment]{SummarizedExperiment-class}}.
#'
#' @param ses A (named) list of
#' \code{\link[SummarizedExperiment]{SummarizedExperiment-class}}.
#' @param genes A vector of genes/row.names to plot.
#' @param do.scale Logical; whether to scale rows in each SE (default TRUE).
#' @param uniqueScale Logical; whether to force the same colorscale for
#' each heatmap.
#' @param assayName The name of the assay to use; if multiple names are given,
#' the first available will be used. Defaults to "logcpm", "lognorm".
#' @param sortBy Names or indexes of `ses` to use for sorting rows (default all)
#' @param only.common Logical; whether to plot only rows common to all SEs
#' (default TRUE).
#' @param cluster_cols Logical; whether to cluster columns (default FALSE).
#' @param cluster_rows Logical; whether to cluster rows (default TRUE if
#' `do.sortRows=FALSE`, FALSE otherwise).
#' @param toporder Optional verctor of categories on which to supra-order when
#' sorting rows, or name of a `rowData` column to use for this purpose.
#' @param hmcols Colors for the heatmap.
#' @param breaks Breaks for the heatmap colors. Alternatively, symmetrical
#' breaks can be generated automatically by setting `breaks` to a numerical
#' value between 0 and 1. The value is passed as the `split.prop` argument to
#' the \code{\link{getBreaks}} function, and indicates the proportion of the
#' points to map to a linear scale, while the more extreme values will be
#' plotted on a quantile scale. `breaks=FALSE` will disable symmetrical scale
#' and quantile capping, while retaining automatic breaks. `breaks=1` will
#' produce a symmetrical scale without quantile capping.
#' @param gaps_at Columns of `colData` to use to establish gaps between columns.
#' @param gaps_row A named vector according to which rows will be split.
#' @param anno_rows Columns of `rowData` to use for annotation.
#' @param anno_columns Columns of `colData` to use for annotation.
#' @param name The title of the heatmap key.
#' @param anno_colors List of colors to use for annotation.
#' @param show_rownames Whether to show row names (default TRUE if 50 rows or
#' less).
#' @param show_colnames Whether to show column names (default FALSE).
#' @param ... Any other parameter passed to each call of
#' \code{\link[ComplexHeatmap]{Heatmap}}.
#'
#' @return A Heatmap list.
#'
#' @importFrom circlize colorRamp2
#' @importFrom methods is
#' @import SummarizedExperiment
#' @import ComplexHeatmap
#' @export
crossHm <- function(ses, genes, do.scale=TRUE, uniqueScale=FALSE,
                    assayName=.getDef("assayName"), sortBy=seq_along(ses),
                    only.common=TRUE, cluster_cols=FALSE,
                    cluster_rows=is.null(sortBy), toporder=NULL, hmcols=NULL,
                    breaks=.getDef("breaks"), gaps_at=.getDef("gaps_at"),
                    gaps_row=NULL, anno_rows=.getDef("anno_rows"),
                    anno_columns=.getDef("anno_columns"), name=NULL,
                    anno_colors=list(), show_rownames=NULL,
                    show_colnames=FALSE, ... ){

    if(is(ses,"SummarizedExperiment")) ses <- list(ses)
    if(is.null(names(ses))) names(ses) <- paste("SE", seq_along(ses))

    tt <- table(unlist(lapply(ses,row.names)))
    if(only.common){
        genes <- intersect(genes, names(tt)[which(tt==length(ses))])
    }else{
        genes <- intersect(genes, names(tt))
    }
    if(length(genes)==0)
        stop("There appears to be not feature in common across `ses`")

    if(!is.null(toporder)){
        if(is.null(names(toporder)))
            stop("`toporder` should be a vector named by feature")
        if(!all(genes %in% names(toporder)))
            warning("Some features are missing from `toporder`")
        toporder <- toporder[genes]
        names(toporder) <- genes
    }

    dats <- lapply( ses, .prepData, genes=genes, assayName=assayName,
                    do.scale=do.scale && !uniqueScale, includeMissing=TRUE )
    if(do.scale && uniqueScale){
        x <- do.call(cbind, dats)
        x <- t(scale(t(x)))
        dl <- lapply(dats, FUN=function(x) seq_len(ncol(x)))
        dl2 <- c(0,cumsum(sapply(dl,length)[-length(dl)]))
        dats <- lapply( seq_along(dl), FUN=function(i)
            x[,dl[[i]]+dl2[[i]],drop=FALSE] )
        names(dats) <- names(dl)
    }else{
        x <- NULL
    }

    if(uniqueScale){
        if(is.null(x)) x <- do.call(cbind, dats)
        if(is.null(breaks) && do.scale) breaks <- 0.995
        cscale <- .prepScale(x, hmcols=.getHMcols(hmcols), breaks=breaks)
        breaks <- cscale$breaks
        hmcols <- cscale$hmcols
    }

    if(!is.null(sortBy) && length(sortBy)>0){
        xs <- dats
        if(do.scale && !uniqueScale)
            xs <- lapply(xs, FUN=function(x){ t(scale(t(x))) })
        xs <- do.call(cbind, xs[sortBy])
        genes <- row.names(sortRows(xs,toporder=toporder))
        dats <- lapply(dats, FUN=function(x) x[genes,,drop=FALSE])
    }

    CDs <- lapply(ses, ac=anno_columns, FUN=function(x,ac){
        ac <- intersect(colnames(colData(x)),ac)
        as.data.frame(colData(x)[,ac,drop=FALSE])
    })
    anno_columns <- colnames(CDs[[1]])

    ses <- lapply(seq_along(ses), FUN=function(i){
        RD <- rowData(ses[[i]])[genes,,drop=FALSE]
        row.names(RD) <- genes
        SummarizedExperiment( list(a=dats[[i]]), colData=CDs[[i]], rowData=RD,
                              metadata=metadata(ses[[i]]) )
    })
    names(ses) <- names(dats)
    anno_colors <- .getAnnoCols(ses[[1]], anno_colors, do.assign=TRUE)

    if(is.null(show_rownames)) show_rownames <- length(genes)<50

    hlp <- list()
    if(uniqueScale){
        if(!is.null(assayName) && length(assayName)==1 &&
           !is.numeric(assayName)){
            hlp$title <- ifelse(do.scale, paste0("scaled\n",assayName), assayName)
        }else{
            hlp$title <- ifelse(do.scale, "z-scores", "")
        }
    }

    htlist <- sapply(seq_along(ses), FUN=function(i){
        sechm(ses[[i]], do.scale=(do.scale && !uniqueScale),
              assayName="a", name=names(ses)[i], toporder=toporder,
              hmcols=hmcols, breaks=breaks, anno_rows=anno_rows,
              anno_columns=anno_columns, anno_colors=anno_colors,
              cluster_rows=FALSE, cluster_cols=cluster_cols, sortRowsOn=NULL,
              show_rownames=(show_rownames && i==length(ses)),
              show_colnames=show_colnames, isMult=i!=length(ses),
              show_heatmap_legend=(!uniqueScale || i==length(ses)),
              heatmap_legend_param=hlp, column_title=names(ses)[i], ...)
    })

    ht <- NULL
    for(f in htlist) ht <- ht + f
    ht
}
