#' meltSE
#'
#' Melts a SE object into a ggplot-compatible data.frame
#'
#' @param x An object of class `SummarizedExperiment`
#' @param genes A vector of genes to include
#' @param assayName The name of the assay to use. If NULL, the first one will be used.
#' @param colDat.columns The data columns to include (defaults includes all)
#' @param value.name The name of the value column.
#'
#' @return A data.frame.
#' @export
meltSE <- function(x, genes=NULL, assayName=NULL, colDat.columns=NULL, value.name=NULL){
  if(is.null(genes)) genes <- row.names(x)
  if(is.null(colDat.columns)) colDat.columns <- colnames(colData(x))
  if(is.null(assayName) && !is.null(assayNames(x)) ) assayName <- assayNames(x)[[1]]
  if(is.null(assayName) || assayName==""){
    a <- assay(x)
    if(is.null(value.name)) value.name <- "value"
  }else{
    a <- assays(x)[[assayName]]
    if(is.null(value.name)) value.name <- assayName
  }
  a <- a[genes,,drop=F]
  df <- data.frame( feature=rep(row.names(a),ncol(a)),
                    sample=rep(colnames(a), each=nrow(a)) )
  for(f in colDat.columns) df[[f]] <- rep(colData(x)[[f]],each=nrow(a))
  df[[value.name]] <- as.numeric(a)
  df
}
