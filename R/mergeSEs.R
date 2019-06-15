#' mergeSEs
#'
#' Merges a list of SummarizedExperiments.
#'
#' @param ll A (named) list of SummarizedExperiments
#' @param what Values to use, either "asis" (input values) or "zscores" (within-dataset row z-scores)
#' @param commonOnly Logical; whether to restrict to rows present in all datasets (default TRUE).
#' @param assayName An optional vector of assayNames to use. The first available will be used, or the first assay if NULL.
#' @param colColumns A character vector specifying `colData` columns to include (if available in at least one of the datasets).
#' If NULL, everything is kept.
#' @param defValues A list specifying the default values for `colColumns` when these are absent.
#'
#' @return An object of class `SummarizedExperiment`
#'
#' @export
mergeSEs <- function(ll, what="zscores", commonOnly=TRUE, assayName=c("logcpm","lognorm"), colColumns=NULL, defValues=list()){
  what <- match.arg(what, c("asis", "zscores"))
  if(!commonOnly && what=="zscores") stop("For z-scores, `commonOnly` must be TRUE.")
  suppressPackageStartupMessages({
    library(SummarizedExperiment)
    library(data.table)
  })
  dats <- lapply(ll, an=assayName, FUN=function(x,an){
    if(is.null(an) || length(an)==0) return(assay(x))
    an <- intersect(an,assayNames(x))
    if(length(an)==0){
      an <- assayNames(x)[1]
      warning(paste0("Specified assay(s) not found in some datasets! Using use the first assay (",an,")."))
    }else{
      an <- an[1]
    }
    assays(x)[[an]]
  })
  tt <- table(unlist(lapply(dats,row.names)))
  if(commonOnly){
    g <- names(tt)[which(tt==length(ll))]
  }else{
    g <- names(tt)
  }
  dats <- lapply(dats, g=g, FUN=function(x,g){
    as.matrix(as.data.frame(x[g,,drop=F]))
  })
  if(what=='zscores'){
    dats <- lapply(dats, FUN=function(y){
      t(apply(y,1,FUN=function(x){
        x2 <- try(as.numeric(scale(as.numeric(x))),silent=T)
        if(is(x2,"try-error")) return(rep(0,length(x)))
        x2
      }))
    })
  }
  dat <- do.call(cbind, dats)

  cd <- lapply(ll, cc=colColumns, FUN=function(x,cc){
    x <- as.data.frame(colData(x))
    if(!is.null(cc)) x <- x[,intersect(cc,colnames(x)),drop=F]
    x
  })
  cd2 <- as.data.frame(data.table::rbindlist(cd, fill=T, idcol="Dataset"))
  for(i in names(defValues)) cd2[[i]][which(is.na(cd2[[i]]))] <- defValues[[i]]

  se <- SummarizedExperiment( dat, colData=cd2)
  nn <- do.call(c, lapply(ll, colnames))
  if(!is.null(names(ll))) nn <- paste(rep(names(ll), sapply(ll,ncol)), nn, sep=".")
  colnames(se) <- nn
  se
}
