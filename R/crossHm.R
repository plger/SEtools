#' crossHm
#'
#' Plot a multi-panel heatmap from a list of
#' \code{\link[SummarizedExperiment]{SummarizedExperiment-class}}
#'
#' @param ses A (named) list of
#' \code{\link[SummarizedExperiment]{SummarizedExperiment-class}}.
#' @param genes A vector of genes/row.names to plot.
#' @param do.scale Logical; whether to scale rows in each SE (default TRUE).
#' @param uniqueColorScale Logical; whether to force the same colorscale for
#' each heatmap (default FALSE; applicable only when what!='asis').
#' @param assayName The name of the assay to use; if multiple names are given,
#' the first available will be used. Defaults to "logcpm", "lognorm".
#' @param do.sortRows Logical; whether to sort rows according to MDS
#' (default TRUE).
#' @param only.common Logical; whether to plot only rows common to all SEs
#' (default TRUE).
#' @param anno_columns A vector of colData columns to use (if available) for
#' annotation, default "Condition" and "TimePoint".
#' @param spreadAnnotation Logical; whether to spread annotation to all SEs
#' (not yet supported)
#' @param hmcols An optional color palette, such as produced by the
#' `colorRampPalette` function.
#' @param cluster_columns Logical; whether to cluster columns (default FALSE).
#' @param cluster_rows Logical; whether to cluster rows (default TRUE if
#' `do.sortRows=FALSE`, FALSE otherwise).
#' @param show_rownames Logical; whether to show row names (default TRUE if
#' `length(genes)<80`, FALSE otherwise).
#' @param show_colnames Logical; whether to show column names (default FALSE)
#' @param anno_colors A vector of color for annotations.
#' @param ... Any other parameter passed to each call of
#' \code{\link[ComplexHeatmap]{Heatmap}}.
#'
#' @return A \code{\link[ComplexHeatmap]{Heatmap-class}} object.
#'
#' @examples
#' data("SE", package="SEtools")
#' crossHm(list(se1=SE[,1:10], se2=SE[,11:20]), head(row.names(SE)))
#'
#' @importFrom circlize colorRamp2
#' @importFrom methods is
#' @import SummarizedExperiment
#' @import ComplexHeatmap
#' @export
crossHm <- function( ses, genes, do.scale=TRUE, uniqueColorScale=FALSE,
                     assayName=.getDef("assayName"), do.sortRows=TRUE,
                     only.common=TRUE, anno_columns=.getDef("anno_columns"),
                     spreadAnnotation=FALSE, hmcols=NULL,
                     cluster_columns=FALSE, cluster_rows=!do.sortRows,
                     show_rownames=ifelse(length(genes)<80,"once",FALSE),
                     show_colnames=FALSE, anno_colors=.getDef("anno_colors"), ...){
  if(is(ses,"SummarizedExperiment")) ses <- list(ses)
  if(is.null(anno_colors)) anno_colors <- list()
  if(only.common){
    tt <- table(unlist(lapply(ses,row.names)))
    genes <- intersect(genes, names(tt)[which(tt==length(ses))])
  }else{
    genes <- intersect(genes, row.names(ses[[1]]))
  }

  if(spreadAnnotation){
    ac <- unlist(lapply(ses, ac=anno_columns, FUN=function(x,ac){
      ac <- intersect(colnames(colData(x)),ac)
      if(length(ac)==0) return(list())
      names(ac) <- ac
      lapply(ac, cd=colData(x), FUN=function(x,cd){ levels(factor(cd[[x]])) })
    }), recursive=FALSE)
    # TO DO
    warning("`spreadAnnotation` not yet implemented.")
  }


  dats <- lapply(names(ses), an=assayName, g=genes, FUN=function(x, g, an){
    dat <- .chooseAssay(ses[[x]], an)
    dat <- as.data.frame(dat)[g,,drop=FALSE]
    if(do.scale){
      dat <- t(apply(dat,1,FUN=function(x){
        tryCatch(x <- scale(x), error=function(e){ })
      }))
    }
    row.names(dat) <- genes
    as.matrix(dat)
  })

  hmcols <- .getHMcols(hmcols)
  if(any(dats[[1]]<0) && uniqueColorScale && !is.function(hmcols)){
    hmcols <- colorRamp2(.getBreaks(dats, length(hmcols)+1), hmcols)
  }

  if(do.sortRows){
    tmp <- do.call(cbind, dats)
    genes <- try(row.names(sortRows(tmp)), silent=TRUE)
    if(is(genes, "try-error")) genes <- row.names(sortRows(tmp, na.rm=TRUE))
    dats <- lapply(dats,g=genes,FUN=function(x,g){ x[g,,drop=FALSE]})
  }

  htlist <- NULL
  for(i in 1:length(ses)){
    cd <- as.data.frame(colData(ses[[i]]))
    cd <- cd[,intersect(colnames(cd), anno_columns),drop=FALSE]
    if(ncol(cd)==0){
      an <- NULL
    }else{
      anc <- anno_colors[which(names(anno_colors) %in% colnames(cd))]
      shown <- (i==length(ses))
      if(length(anc)==0){
          an <- HeatmapAnnotation(df=cd, show_annotation_name=shown)
      }else{
          an <- HeatmapAnnotation(df=cd, col=anc, show_annotation_name=shown)
      }
    }
    srn <- show_rownames
    if(srn=="once") srn <- i==length(ses)
    if(!is.function(hmcols)){
        tcols <- colorRamp2(.getBreaks(dats[[i]]),hmcols)
    }else{
        tcols <- hmcols
    }
    legtitle <- ifelse(uniqueColorScale | length(ses)==1, assayName[1], names(ses)[[i]])
    htlist <- htlist + Heatmap(dats[[i]], col=tcols, na_col="white",
                               name=names(ses)[[i]],
                               column_title=names(ses)[[i]], top_annotation=an,
                               cluster_rows=cluster_rows,
                               cluster_columns=cluster_columns,
                               show_row_names=srn,
                               show_heatmap_legend=!(i>1 & uniqueColorScale),
                               heatmap_legend_param=list(title=legtitle),
                               show_column_names=show_colnames,
                               ...)
  }
  htlist
}
