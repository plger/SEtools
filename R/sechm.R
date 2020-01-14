#' Heatmap wrappers for
#' \code{\link[SummarizedExperiment]{SummarizedExperiment-class}}.
#' @name SE-heatmap
#' @rdname SE-heatmap
#' @aliases sehm sechm
#'
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
#' @param gaps_row Passed to \code{\link[pheatmap]{pheatmap}}; if missing, will
#' be set automatically according to toporder.
#' @param anno_rows Columns of `rowData` to use for annotation.
#' @param anno_columns Columns of `colData` to use for annotation.
#' @param anno_colors List of colors to use for annotation.
#' @param show_rownames Whether to show row names (default TRUE if 50 rows or
#' less).
#' @param show_colnames Whether to show column names (default FALSE).
#' @param ... Further arguments passed to `pheatmap` (`sehm`) or `Heatmap`
#' (`sechm`).
#'
#' @return A heatmap (see \code{\link[pheatmap]{pheatmap}}), or, for `sechm`,
#' a \code{\link[ComplexHeatmap]{Heatmap-class}}.
#'
#' @examples
#' data("SE", package="SEtools")
#' sehm(SE, do.scale=TRUE)

#' @rdname SE-heatmap
#' @param isMult Logical; used to silence labels when plotting mulitple heatmaps
#' @param show_heatmap_legend Logical; whether to show heatmap legend
#'
#' @importFrom circlize colorRamp2
#' @importFrom methods is
#' @import SummarizedExperiment
#' @import ComplexHeatmap
#' @export
sechm <- function(se, genes=NULL, do.scale=FALSE, assayName=.getDef("assayName"),
                  sortRowsOn=seq_len(ncol(se)), cluster_cols=FALSE,
                  cluster_rows=is.null(sortRowsOn), toporder=NULL, hmcols=NULL,
                  breaks=.getDef("breaks"), gaps_at=.getDef("gaps_at"),
                  gaps_row=NULL, anno_rows=.getDef("anno_rows"),
                  anno_columns=.getDef("anno_columns"), name=NULL,
                  anno_colors=list(), show_rownames=NULL, show_colnames=FALSE,
                  isMult=FALSE, show_heatmap_legend=!isMult,
                  ...){

  assayName <- .chooseAssay(se, assayName, returnName = TRUE)
  if(is.null(name)){
      if(is.numeric(assayName)){
          name <- ifelse(do.scale, "z-scores", "assay")
      }else{
          name <- ifelse(do.scale, paste0("scaled\n",assayName), assayName)
      }
  }

  x <- .prepData(se, genes=genes, do.scale=do.scale, assayName=assayName)

  toporder <- .parseToporder(rowData(se)[row.names(x),], toporder)
  if(!is.null(sortRowsOn) || length(sortRowsOn)>0)
    x <- x[row.names(sortRows(x[,sortRowsOn],toporder=toporder,na.rm=TRUE)),]

  if( is.null(breaks) ){
      if( (!is.null(assayName) && grepl("^log[2]?FC$",assayName)) || do.scale)
          breaks <- 0.995
  }
  hmcols <- .getBaseHMcols(se, hmcols)
  cscale <- .prepScale(x, hmcols=hmcols, breaks=breaks)
  breaks <- cscale$breaks
  hmcols <- circlize::colorRamp2(breaks[-length(breaks)], cscale$hmcols)

  anno_colors <- .getAnnoCols(se, anno_colors)
  anr <- .prepareAnnoDF( rowData(se)[row.names(x),], anno_colors,
                         anno_rows, whichComplex="row" )

  an <- .prepareAnnoDF( colData(se), anno_colors, anno_columns,
                        whichComplex="column", show_legend=!isMult,
                        show_annotation_name=!isMult )

  gaps_col <- gaps_row <- NULL
  if(!is.null(gaps_at) &&
     length(gaps_at <- intersect(gaps_at, colnames(colData(se))))>0){
      gaps_col <- as.data.frame(colData(se))[,gaps_at,drop=FALSE]
  }
  if(!is.null(gaps_row) &&
     length(gaps_row <- intersect(gaps_row, colnames(rowData(se))))>0){
      gaps_row <- as.data.frame(rowData(se))[x,gaps_row,drop=FALSE]
  }

  if(is.null(show_rownames)) show_rownames <- nrow(x)<50
  Heatmap(x, col=hmcols, na_col="white", name=name,
          show_row_names=show_rownames, show_column_names=show_colnames,
          top_annotation=an, left_annotation=anr, row_split=gaps_row,
          column_split=gaps_col, show_heatmap_legend=show_heatmap_legend,
          cluster_rows=cluster_rows, cluster_columns=cluster_cols,
          ...)
}
