#' sehm
#'
#' Heatmap wrapper for \code{\link[SummarizedExperiment]{SummarizedExperiment-class}}.
#'
#' @param se A \code{\link[SummarizedExperiment]{SummarizedExperiment-class}}.
#' @param genes An optional vector of genes (i.e. row names of `se`)
#' @param do.scale Logical; whether to scale rows (default FALSE).
#' @param assayName An optional vector of assayNames to use. The first available
#'  will be used, or the first assay if NULL.
#' @param sortRowsOn Sort rows by MDS polar order using the specified columns
#' (default all)
#' @param cluster_cols Whether to cluster columns (default F)
#' @param cluster_rows Whether to cluster rows; default FALSE if
#' `do.sortRows=TRUE`.
#' @param toporder Optional verctor of categories on which to supra-order when
#' sorting rows.
#' @param hmcols Colors for the heatmap.
#' @param breaks Breaks for the heatmap colors. Alternatively, if
#' `breaks==TRUE`, a symmetrical scale with capped ends will be used
#' (appropriate when plotting log2 foldchanges)
#' @param gaps_at Columns of `colData` to use to establish gaps between columns.
#' @param anno_rows Columns of `rowData` to use for annotation.
#' @param anno_columns Columns of `colData` to use for annotation.
#' @param anno_colors List of colors to use for annotation.
#' @param show_rownames Whether to show row names (default TRUE if 50 rows or
#' less).
#' @param show_colnames Whether to show column names (default FALSE).
#' @param ... Further arguments passed to `pheatmap`.
#'
#' @return A heatmap (see \code{\link[pheatmap]{pheatmap}}).
#'
#' @examples
#' data("SE", package="SEtools")
#' sehm(SE, do.scale=TRUE)
#'
#' @importFrom pheatmap pheatmap
#' @import SummarizedExperiment
#' @export
sehm <- function( se, genes=NULL, do.scale=FALSE, assayName=.getDef("assayName"),
                  sortRowsOn=1:ncol(se), cluster_cols=FALSE,
                  cluster_rows=is.null(sortRowsOn), toporder=NULL, hmcols=NULL,
                  breaks=.getDef("breaks"), gaps_at=.getDef("gaps_at"),
                  anno_rows=.getDef("anno_rows"),
                  anno_columns=.getDef("anno_columns"),
                  anno_colors=.getDef("anno_colors"), show_rownames=NULL,
                  show_colnames=FALSE, ...){

  x <- as.matrix(.chooseAssay(se, assayName))

  if(!is.null(genes)) x <- x[intersect(genes,row.names(x)),]
  if(do.scale){
    x <- x[apply(x,1,FUN=sd)>0,]
    x <- t(scale(t(x)))
  }
  if(!is.null(sortRowsOn)){
    if(!is.null(toporder)){
      if(!is.null(names(toporder))){
        toporder <- toporder[row.names(x)]
      }
    }
    x2 <- sortRows(x[,sortRowsOn], toporder=toporder, na.rm=TRUE)
    x <- x[row.names(x2),]
    rm(x2)
  }

  hmcols <- .getHMcols(hmcols)
  if(!is.null(breaks) && length(breaks)==1 && breaks==TRUE)
    breaks <- .getBreaks(x, length(hmcols))

  anr <- as.data.frame(rowData(se))
  anr <- anr[,intersect(colnames(anr), anno_rows),drop=FALSE]
  if(ncol(anr)==0) anr <- NULL
  an <- as.data.frame(colData(se))
  an <- an[,intersect(colnames(an), anno_columns),drop=FALSE]
  if(ncol(an)==0){
    an <- NULL
  }else{
    for(i in colnames(an)){
      if(is.logical(an[[i]])){
        an[[i]] <- factor(as.character(an[[i]]),levels=c("FALSE","TRUE"))
        if(!(i %in% names(anno_colors))){
          anno_colors[[i]] <- c("FALSE"="white", "TRUE"="darkblue")
        }
      }else{
        if(is.factor(an[[i]])) an[[i]] <- droplevels(an[[i]])
      }
    }
  }


  if(!is.null(gaps_at)){
    gaps_at <- intersect(gaps_at, colnames(colData(se)))
    ga <- apply( as.data.frame(colData(se))[,gaps_at,drop=FALSE], 1,
                 collapse=" ", FUN=paste)
    ga <- factor(ga, levels=unique(ga))
    o <- order(ga)
    x <- x[,o]
    an <- an[o,,drop=FALSE]
    ga <- ga[o]
    gaps <- (which(!duplicated(ga))-1)[-1]
  }else{
    gaps <- NULL
  }

  if(is.null(show_rownames)) show_rownames <- nrow(x) <= 50
  pheatmap(x, color=hmcols, border_color=NA, gaps_col=gaps, breaks=breaks,
           cluster_cols=cluster_cols, cluster_rows=cluster_rows,
           annotation_col=an, annotation_row=anr, annotation_colors=anno_colors,
           show_rownames=show_rownames, show_colnames=show_colnames, ...)
}
