#' sehm
#'
#' Deprecated pheatmap wrapper for
#' \code{\link[SummarizedExperiment]{SummarizedExperiment-class}}.
#' **This function has been replaced by the \code{\link[sechm]{sechm}} function
#' from the `sechm` package and is retained here solely for backward 
#' compatibility.**
#'
#' @param se A \code{\link[SummarizedExperiment]{SummarizedExperiment-class}}.
#' @param genes An optional vector of genes (i.e. row names of `se`)
#' @param do.scale Logical; whether to scale rows (default FALSE).
#' @param assayName An optional vector of assayNames to use. The first available
#'  will be used, or the first assay if NULL.
#' @param cluster_cols Whether to cluster columns (default F)
#' @param cluster_rows Whether to cluster rows; default FALSE if
#' `do.sortRows=TRUE`.
#' @param sortRowsOn Sort rows by MDS polar order using the specified columns
#' (default all)
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
#' @param gaps_row Passed to the heatmap function; if missing, will
#' be set automatically according to toporder.
#' @param anno_rows Columns of `rowData` to use for left annotation.
#' @param anno_columns Columns of `colData` to use for top annotation.
#' @param anno_colors List of colors to use for annotation.
#' @param show_rownames Whether to show row names (default TRUE if less than
#' 50 rows to plot).
#' @param show_colnames Whether to show column names (default FALSE).
#' @param ... Further arguments passed to `pheatmap`
#'
#' @return A heatmap.
#'
#' @importFrom circlize colorRamp2
#' @importFrom methods is
#' @import SummarizedExperiment pheatmap sechm
#' @export
sehm <- function( se, genes, do.scale=FALSE, assayName=.getDef("assayName"),
                  sortRowsOn=seq_len(ncol(se)), cluster_cols=FALSE,
                  cluster_rows=is.null(sortRowsOn), toporder=NULL, hmcols=NULL,
                  breaks=.getDef("breaks"), gaps_at=.getDef("gaps_at"),
                  gaps_row=NULL, anno_rows=.getDef("anno_rows"),
                  anno_columns=.getDef("anno_columns"),
                  anno_colors=NULL, show_rownames=NULL,
                  show_colnames=FALSE, ...){
  .Deprecated(msg="'sehm' is deprecated; please use sechm::sechm instead")

  if(is.null(anno_colors)) anno_colors <- sechm:::.getAnnoCols(se)
  assayName <- sechm:::.chooseAssay(se, intersect(assayName,names(assays(se))),
                                    returnName=TRUE)
  x <- sechm:::.prepData(se, genes=genes, do.scale=do.scale, assayName=assayName)
  
  toporder <- sechm:::.parseToporder(rowData(se)[row.names(x),,drop=FALSE], toporder)
  if(!is.null(sortRowsOn) && length(sortRowsOn)>0 && nrow(x)>2){
    x <- x[row.names(sechm::sortRows(x[,sortRowsOn,drop=FALSE],
                                     toporder=toporder,na.rm=TRUE)),]
  }
  if( is.null(breaks) ){
    if( (!is.null(assayName) && grepl("^log[2]?FC$",assayName)) || do.scale)
      breaks <- 0.995
  }
  hmcols <- sechm:::.getBaseHMcols(se, hmcols)
  cscale <- sechm:::.prepScale(x, hmcols=hmcols, breaks=breaks)
  breaks <- cscale$breaks
  hmcols <- cscale$hmcols
  
  anr <- sechm:::.prepareAnnoDF( rowData(se)[row.names(x),,drop=FALSE], 
                                 anno_colors, anno_rows )
  anno_colors <- anr$anno_colors
  anr <- anr$an
  
  an <- sechm:::.prepareAnnoDF( colData(se), anno_colors, anno_columns )
  anno_colors <- an$anno_colors
  an <- an$an
  
  gaps_col <- sechm:::.getGaps(gaps_at, colData(se), silent=TRUE)
  if(!is.null(gaps_row) && !is.logical(gaps_row))
    gaps_row <- .getGaps(gaps_row, rowData(se)[row.names(x),,drop=FALSE])
  
  if(!is.null(gaps_col)){
    ga <- apply( gaps_col, 1, collapse=" ", FUN=paste)
    ga <- factor(ga, levels=unique(ga))
    o <- order(ga)
    x <- x[,o]
    an <- an[o,,drop=FALSE]
    ga <- ga[o]
    gaps_col <- (which(!duplicated(ga))-1)[-1]
  }
  
  if(!is.null(toporder) && is.null(gaps_row)){
    toporder <- toporder[row.names(x)]
    gaps_row <- (which(!duplicated(toporder))-1)[-1]
  }else if(!is.null(gaps_row) && !all(gaps_row==FALSE)){
    ga <- apply( gaps_row, 1, collapse=" ", FUN=paste)
    ga <- factor(ga, levels=unique(ga))
    o <- order(ga)
    x <- x[o,]
    anr <- anr[o,,drop=FALSE]
    ga <- ga[o]
    gaps_row <- (which(!duplicated(ga))-1)[-1]
  }
  
  if(is.null(show_rownames)) show_rownames <- nrow(x) <= 50
  if(nrow(x)<=2) cluster_rows <- FALSE
  pheatmap::pheatmap(x, color=hmcols, border_color=NA, gaps_col=gaps_col, breaks=breaks,
                     gaps_row=gaps_row, cluster_cols=cluster_cols, cluster_rows=cluster_rows,
                     annotation_col=an, annotation_row=anr, annotation_colors=anno_colors,
                     show_rownames=show_rownames, show_colnames=show_colnames, ...)
}

.getDef <- function(x, ...) sechm:::.options$get(x)
