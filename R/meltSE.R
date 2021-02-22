#' meltSE
#'
#' Melts a SE object into a \code{\link[ggplot2]{ggplot}}-ready long data.frame.
#'
#' @param x An object of class
#' \code{\link[SummarizedExperiment]{SummarizedExperiment-class}}
#' @param genes A vector of genes to include. Use `genes=NULL` to include all.
#' @param assayName The name(s) of the assay(s) to use. If NULL and the assays are named,
#' all of them will be included (if they are not named, the first one will be used).
#' @param colDat.columns The colData columns to include (defaults includes all).
#' Use `colDat.columns=NA` in order not to include any.
#' @param rowDat.columns The rowData columns to include (none included by default). Use
#' `rowData=NULL` to include all.
#'
#' @return A data.frame.
#'
#' @examples
#' data("SE", package="SEtools")
#' head(meltSE(SE,"Fos"))
#'
#' @import SummarizedExperiment
#' @export
meltSE <- function(x, genes, assayName=NULL, colDat.columns=NULL,
                   rowDat.columns=NA){
  genes <- intersect(genes, row.names(x))
  if(is.null(colDat.columns)) colDat.columns <- colnames(colData(x))
  if(all(is.na(colDat.columns))) colDat.columns <- c()
  colDat.columns <- intersect(colDat.columns, colnames(colData(x)))
  if(is.null(rowDat.columns)) rowDat.columns <- colnames(rowData(x))
  if(all(is.na(rowDat.columns))) rowDat.columns <- c()
  rowDat.columns <- intersect(rowDat.columns, colnames(rowData(x)))
  if(is.null(assayName) && !is.null(assayNames(x)) ) assayName <- assayNames(x)
  if(is.null(assayName)){
    a <- list(value=assay(x))
  }else{
    a <- assays(x)[assayName]
  }
  if(is.numeric(assayName)) names(a) <- paste0("assay", assayName)
  a <- lapply(a, FUN=function(x) x[genes,,drop=FALSE])
  df <- data.frame( feature=rep(genes,ncol(x)),
                    sample=rep(colnames(x), each=length(genes)) )
  for(f in colDat.columns) df[[f]] <- rep(colData(x)[[f]],each=length(genes))
  for(f in rowDat.columns) df[[f]] <- rep(rowData(x)[genes,f], ncol(x))
  for(v in names(a)) df[[v]] <- as.numeric(a[[v]])
  df
}



#' castSE
#'
#' Casts a data.frame as a \code{\link[SummarizedExperiment]{SummarizedExperiment-class}}
#'
#' @param x A data.frame
#' @param rowNames Column of `x` containing the row.names (if omitted, will build from
#' `rowData`)
#' @param colNames Column of `x` containing the column names (if omitted, will build from
#' `colData`)
#' @param assayNames Columns of `x` to turn into assays
#' @param colData Columns of `x` to use as colData
#' @param rowData Columns of `x` to use as rowData
#' @param sparse Local, whether to keep the assays sparse.
#'
#' @return A \code{\link[SummarizedExperiment]{SummarizedExperiment-class}}
#' @export
#' @import Matrix
#' @importFrom stats setNames
#'
#' @examples
#' d <- data.frame(transcript=rep(LETTERS[1:10],each=2), gene=rep(LETTERS[1:5],each=4),
#'                 count=rpois(20, 10), sample=letters[1:2])
#' head(d)
#' castSE(d, rowData=c("transcript","gene"), colNames="sample")
castSE <- function(x, rowNames=NULL, colNames=NULL, assayNames=NULL, colData=NULL, rowData=NULL, sparse=FALSE){
    if(!is.data.frame(x) && !is(x,"DFrame")) stop("`x` should be a data frame.")
    rowsDef <- !(is.null(rowNames) && is.null(rowData))
    colsDef <- !(is.null(colNames) && is.null(colData))
    if(!rowsDef & !colsDef) stop("Insufficient information provided.")
    if(is.null(assayNames)){
      assayNames <- setdiff(colnames(x), c(colNames, rowNames, colData, rowData))
    }
    if(colsDef){
      if(is.null(colNames)){
          colNames <- as.factor(do.call(paste0, c(list(sep="."),lapply(colData, FUN=function(y) x[[y]]))))
      }else{
          colNames <- as.factor(x[[colNames]])
      }
    }
    if(rowsDef){
        if(is.null(rowNames)){
            rowNames <- as.factor(do.call(paste, c(list(sep="."),lapply(rowData, FUN=function(y) x[[y]]))))
        }else{
            rowNames <- as.factor(x[[rowNames]])
        }
    }
    rowNames <- droplevels(rowNames)
    colNames <- droplevels(colNames)
    ass <- lapply(setNames(assayNames, assayNames), FUN=function(y){
        sparseMatrix(i=as.integer(rowNames), j=as.integer(colNames), x=x[[y]],
                     dim=c(length(levels(rowNames)), length(levels(colNames))),
                    dimnames=list(levels(rowNames), levels(colNames)))
    })
    rowData <- x[!duplicated(rowNames),rowData,drop=FALSE]
    row.names(rowData) <- rowNames[!duplicated(rowNames)]
    colData <- x[!duplicated(colNames),colData,drop=FALSE]
    row.names(colData) <- colNames[!duplicated(colNames)]
    if(!sparse) ass <- lapply(ass, as.matrix)
    SummarizedExperiment(ass, colData=colData[colnames(ass[[1]]),,drop=FALSE],
                         rowData=rowData[row.names(ass[[1]]),,drop=FALSE])
}
