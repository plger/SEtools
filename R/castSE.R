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
