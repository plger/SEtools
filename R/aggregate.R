.aggRowDat <- function(x){
    if(is.logical(x)) return(sum(x,na.rm=TRUE)/length(x))
    if(is.numeric(x)) return(median(x,na.rm=TRUE))
    if(is.factor(x)) x <- as.character(x)
    x <- x[!is.na(x)]
    if(length(x)==0) return(NA_character_)
    paste(sort(unique(x)), collapse=";")
}

#' aggSE
#'
#' Aggregates the rows of a `SummarizedExperiment`.
#'
#' @param x An object of class `SummarizedExperiment`
#' @param by Vector by which to aggregate, or column of `rowData(x)`
#' @param assayFun Function by which to aggregate, or a list of such functions
#' (or vector of function names) of the same length as there are assays. If NULL
#' will attempt to use an appropriate function (and notify the functions used),
#' typically the mean.
#' @param rowDatFun Function by which to aggregate the rowData; by default,
#' logical are transformed into a proportion, numerics are aggregated by median,
#' and unique factors/characters are pasted together.
#'
#' @return An object of class `SummarizedExperiment`
#' @export
#'
#' @import SummarizedExperiment
#' @examples
#' data("SE", package="SEtools")
#' # arbitrary IDs for example aggregation:
#' rowData(SE)$otherID <- rep(LETTERS[1:10],each=10)
#' SE <- aggSE(SE, "otherID")
aggSE <- function(x, by, assayFun=NULL, rowDatFun=.aggRowDat){
    if(!is(x,"SummarizedExperiment"))
        stop("`x` should be a `SummarizedExperiment`.")
    if(is.character(by) && length(by)==1){
        if(!(by %in% colnames(rowData(x))))
            stop(paste0("`",by,"` not found in rowData!"))
        by <- rowData(x)[[by]]
    }
    if(is.null(assayFun)){
        if(is.null(assayNames(x))){
            message("`assayFun` undefined and no assayNames in object: will aggregate using means.")
            assayFun <- mean
        }else{
            sumassays <- c("counts","cpm","tpm","rpkm","fpkm")
            assayFun <- sapply(assayNames(x), FUN=function(x){
                if(tolower(x) %in% sumassays) return("sum")
                if(tolower(x) %in% paste0("log",sumassays)) return("expsum")
                return("mean")
            })
            message("Aggregation methods for each assay:\n",
                    paste(paste0(assayNames(x), ": ", assayFun), collapse="; "))
        }
    }else{
        if(length(assayFun)==1){
            assayFun <- lapply(assays(x), FUN=function(x) assayFun)
        }else{
            if(length(assayFun) != length(assays(x)))
                stop("length(assayFun) != length(assays(x))")
        }
    }
    if(is.null(assayNames(x))) assayNames(x) <- paste0("A", seq_along(assays(x)))
    if(is.null(names(assayFun))) names(assayFun) <- assayNames(x)
    assayFun <- assayFun[assayNames(x)]
    a <- lapply(assayNames(x), FUN=function(y){
        agf <- assayFun[[y]]
        if(agf=="expsum") agf <- function(x) log(sum(exp(x)))
        y <- assays(x)[[y]]
        y <- aggregate(y, by=list(by), FUN=agf)
        row.names(y) <- y[,1]
        as.matrix(y[,-1])
    })
    names(a) <- assayNames(x)
    if(is.null(rowDatFun)){
        return(SummarizedExperiment(a, colData=colData(x), metadata=x@metadata))
    }
    rd <- as.data.frame(rowData(x))
    for(f in colnames(rd)) if(is.factor(rd[[f]])) rd[[f]] <- as.character(rd[[f]])
    rd <- aggregate(rowData(x), list(by), FUN=rowDatFun)
    row.names(rd) <- rd[,1]
    rd <- rd[,-1]
    SummarizedExperiment(a, colData=colData(x), rowData=rd, metadata=x@metadata)
}
