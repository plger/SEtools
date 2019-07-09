#' crossHm
#'
#' Plot a multi-panel heatmap from a list of `SummarizedExperiment`
#'
#' @param ses A (named) list of `SummarizedExperiment` (SEs) objects.
#' @param genes A vector of genes/row.names to plot.
#' @param what What to plot; either "zscores" (within-SE z-scores, default), "asis" (the input data), or "log2FC".
#' If "log2FC".
#' @param uniqueColorScale Logical; whether to force the same colorscale for each heatmap (default FALSE; applicable only when what!='asis').
#' @param ctrlCondition The condition to compare to if `what="log2FC"` (default "Homecage").
#' @param assayName The name of the assay to use; if multiple names are given, the first available will be used.
#' Defaults to "logcpm", "lognorm".
#' @param do.sortRows Logical; whether to sort rows according to MDS (default TRUE).
#' @param only.common Logical; whether to plot only rows common to all SEs (default TRUE).
#' @param acolumns A vector of colData columns to use (if available) for annotation, default "Condition" and "TimePoint".
#' @param spreadAnnotation Logical; whether to spread annotation to all SEs (not yet supported)
#' @param cluster_columns Logical; whether to cluster columns (default FALSE).
#' @param cluster_rows Logical; whether to cluster rows (default TRUE if `do.sortRows=FALSE`, FALSE otherwise).
#' @param show_row_names Logical; whether to show row names (default TRUE if `length(genes)<80`, FALSE otherwise).
#' @param show_column_names Logical; whether to show column names (default FALSE)
#' @param acolors A vector of color for annotations.
#' @param ... Any other parameter passed to each call of `Heatmap`.
#'
#' @importFrom circlize colorRamp2
#' @import SummarizedExperiment
#' @import ComplexHeatmap
#' @export
crossHm <- function( ses, genes, what="zscores", uniqueColorScale=FALSE, ctrlCondition="Homecage", assayName=c("logcpm","lognorm"),
                     do.sortRows=TRUE, only.common=TRUE, acolumns=c("Condition","TimePoint"),
                     spreadAnnotation=FALSE, cluster_columns=FALSE, cluster_rows=!do.sortRows,
                     show_row_names=ifelse(length(genes)<80,"once",FALSE), show_column_names=FALSE, acolors=NULL, ...){
  what <- match.arg(what, c("asis", "zscores", "log2FC"))
  if(is(ses,"SummarizedExperiment")) ses <- list(ses)
  if(is.null(acolors)) acolors <- list()
  if(only.common){
    tt <- table(unlist(lapply(ses,row.names)))
    genes <- intersect(genes, names(tt)[which(tt==length(ses))])
  }else{
    genes <- intersect(genes, row.names(ses[[1]]))
  }

  if(spreadAnnotation){
    ac <- unlist(lapply(ses, ac=acolumns, FUN=function(x,ac){
      ac <- intersect(colnames(colData(x)),ac)
      if(length(ac)==0) return(list())
      names(ac) <- ac
      lapply(ac, cd=colData(x), FUN=function(x,cd){ levels(factor(cd[[x]])) })
    }), recursive=F)
    # TO DO
    warning("`spreadAnnotation` not yet implemented.")
  }


  dats <- lapply(names(ses), an=assayName, g=genes, FUN=function(x, g, an){
    dat <- .chooseAssay(ses[[x]], an)
    dat <- as.data.frame(dat)[g,,drop=FALSE]

    if(what=="zscores"){
      dat <- t(apply(dat,1,FUN=function(x){
        tryCatch(x <- scale(x), error=function(e){ })
      }))
    }
    if(what=="log2FC"){
      w <- which(ses[[i]]$Condition==ctrlCondition)
      dat <- t(apply(dat,1,w=w,FUN=function(x,w){
        cm <- mean(as.numeric(x[w]),na.rm=T)
        if(is.na(cm) || cm==0) return(log2(x+0.1))
        log2(x/cm + 0.1)
      }))
    }
    row.names(dat) <- genes
    dat
  })

  if(what=="asis" || !uniqueColorScale){
    hmcols <- colorRampPalette(c("blue", "black", "yellow"))(29)
  }else{
    vr <- range(unlist(lapply(dats, range)))
    if(max(abs(vr))>3){
      vr <- c(max(abs(vr)), 3, 2.5, 1.5, 1, 0.5)
    }else{
      vr <- c(3, 2.5, 1.5, 1, 0.5)
    }
    vr <- c(-1*vr, 0, rev(vr))
    hmcols <- colorRamp2(vr, colorRampPalette(c("blue", "black", "yellow"))(length(vr)))
  }

  if(do.sortRows){
    tmp <- do.call(cbind, dats)
    genes <- try(row.names(sortRows(tmp)), silent=T)
    if(is(genes, "try-error")) genes <- row.names(sortRows(tmp, na.rm=T))
    dats <- lapply(dats,g=genes,FUN=function(x,g){ x[g,,drop=F]})
  }

  htlist <- NULL
  for(i in 1:length(ses)){
    cd <- as.data.frame(colData(ses[[i]]))
    cd <- cd[,intersect(colnames(cd), acolumns),drop=F]
    if(ncol(cd)==0){
      an <- NULL
    }else{
      an <- HeatmapAnnotation(df=cd, col=acolors[which(names(acolors) %in% colnames(cd))])
    }
    srn <- show_row_names
    if(srn=="once") srn <- i==length(ses)
    htlist <- htlist + Heatmap(dats[[i]], col=hmcols, na_col="white", name=names(ses)[[i]],
                               column_title=names(ses)[[i]], top_annotation=an, cluster_rows=cluster_rows,
                               cluster_columns = cluster_columns, show_row_names=srn,
                               show_column_names=show_column_names, ...)
  }
  htlist
}
