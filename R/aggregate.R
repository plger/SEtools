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
#' @param rowDatFuns A named list providing functions by which to aggregate each
#' rowData columns. If a given column has no specified function, the default
#' will be used, i.e. logical are transformed into a proportion, numerics are
#' aggregated by median, and unique factors/characters are pasted together. Use
#' `rowDataFuns=NULL` to discard rowData.
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
aggSE <- function(x, by, assayFun=NULL, rowDatFuns=list()){
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
    if(is.null(rowDatFuns)){
        return(SummarizedExperiment(a, colData=colData(x), metadata=x@metadata))
    }
    rd <- .aggRowDat(rowData(x), by, agFuns=rowDatFuns)
    SummarizedExperiment(a, colData=colData(x), rowData=rd, metadata=x@metadata)
}

.aggRowDat <- function(rd, by, agFuns=list()){
    if(is.null(agFuns)) agFuns <- list()
    i <- split(seq_len(nrow(rd)), by)
    names(ff) <- ff <- colnames(rd)
    a <- as.data.frame(lapply(ff, FUN=function(y){
        if(y %in% names(agFuns))
            return(sapply(i, FUN=function(x) agFuns[[y]](rd[[y]][x])))
        if(is.logical(rd[[y]]))
            return(vapply(i, FUN.VALUE=vector("numeric",1), FUN=function(x){
                x <- rd[[y]][x]
                sum(x,na.rm=TRUE)/length(x)
            }))
        if(is.numeric(rd[[y]]))
            return(vapply(i, FUN.VALUE=vector(mode(rd[[y]]),1), FUN=function(x){
                x <- rd[[y]][x]
                median(x,na.rm=TRUE)
            }))
        if(is.factor(rd[[y]])) rd[[y]] <- as.character(rd[[y]])
        vapply(i, FUN.VALUE=vector("character",1), FUN=function(x){
            x <- rd[[y]][x]
            x <- x[!is.na(x)]
            if(length(x)==0) return(NA_character_)
            paste(sort(unique(x)), collapse=";")
        })
    }))
    ff <- colnames(rd)[sapply(colnames(rd), FUN=function(x) is.factor(rd[[x]]))]
    for( f in ff ){
        if(all(a[[f]] %in% levels(rd[[f]]))){
            a[[f]] <- factor(a[[f]], levels(rd[[f]]))
        }else{
            a[[f]] <- as.factor(a[[f]])
        }
    }
    a
}
