#' @rdname SE-heatmap
#' @importFrom pheatmap pheatmap
#' @import SummarizedExperiment
#' @importFrom pheatmap pheatmap
#' @export
sehm <- function( se, genes=NULL, do.scale=FALSE, assayName=.getDef("assayName"),
                  sortRowsOn=seq_len(ncol(se)), cluster_cols=FALSE,
                  cluster_rows=is.null(sortRowsOn), toporder=NULL, hmcols=NULL,
                  breaks=.getDef("breaks"), gaps_at=.getDef("gaps_at"),
                  gaps_row=NULL, anno_rows=.getDef("anno_rows"),
                  anno_columns=.getDef("anno_columns"),
                  anno_colors=.getAnnoCols(se), show_rownames=NULL,
                  show_colnames=FALSE, ...){
  ## see sechm.R for a definition of the arguments
  x <- .prepData(se, genes=genes, do.scale=do.scale, assayName=assayName)

  toporder <- .parseToporder(rowData(se)[row.names(x),], toporder)
  if(!is.null(sortRowsOn) && length(sortRowsOn)>0){
      x2 <- sortRows(x[,sortRowsOn,drop=FALSE],toporder=toporder,na.rm=TRUE)
      x <- x[row.names(x2),]
  }
  if( is.null(breaks) ){
      if( (!is.null(assayName) && grepl("^log[2]?FC$",assayName)) || do.scale)
          breaks <- 0.995
  }
  hmcols <- .getBaseHMcols(se, hmcols)
  cscale <- .prepScale(x, hmcols=hmcols, breaks=breaks)
  breaks <- cscale$breaks
  hmcols <- cscale$hmcols

  anr <- .prepareAnnoDF( rowData(se)[row.names(x),], anno_colors, anno_rows )
  anno_colors <- anr$anno_colors
  anr <- anr$an

  an <- .prepareAnnoDF( colData(se), anno_colors, anno_columns )
  anno_colors <- an$anno_colors
  an <- an$an

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
  if(!is.null(toporder) && is.null(gaps_row)){
      toporder <- toporder[row.names(x)]
      gaps_row <- (which(!duplicated(toporder))-1)[-1]
  }else{
      if(length(gaps_row)==1 && is.na(gaps_row)) gaps_row <- NULL
  }

  if(is.null(show_rownames)) show_rownames <- nrow(x) <= 50

  pheatmap(x, color=hmcols, border_color=NA, gaps_col=gaps, gaps_row=gaps_row,
           breaks=breaks, cluster_cols=cluster_cols, cluster_rows=cluster_rows,
           annotation_col=an, annotation_row=anr, annotation_colors=anno_colors,
           show_rownames=show_rownames, show_colnames=show_colnames, ...)
}
