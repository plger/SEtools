#' sortRows
#'
#' @param x A numeric matrix or data.frame.
#' @param z Whether to scale rows for the purpose of calculating order (default FALSE).
#' @param toporder Optional verctor of categories (length=nrow(x)) on which to
#' supra-order  when sorting rows.
#' @param na.rm Wheter to remove missing values and invariant rows (default FALSE).
#' @param method Seriation method; 'MDS_angle' (default) or 'R2E' recommended.
#' @param toporder.meth Whether to perform higher-order sorting 'before'
#' (default) or 'after' the lower-order sorting.
#'
#' @return A reordered matrix or data.frame.
#'
#' @examples
#' # random data
#' m <- matrix( round(rnorm(100,mean=10, sd=2)), nrow=10,
#'              dimnames=list(LETTERS[1:10], letters[11:20]) )
#' m
#' sortRows(m)
#'
#' @importFrom seriation seriate get_order
#' @export
sortRows <- function(x, z=FALSE, toporder=NULL, na.rm=FALSE, method="MDS_angle",
                     toporder.meth="before"){
  toporder.meth <- match.arg(toporder.meth, c("before","after"))
  if(is.numeric(toporder)) toporder <- as.character(toporder)
  if(na.rm){
    w <- which( apply(x, 1, FUN = function(y){ !any(is.na(y)) }) |
                  !(apply(x, 1, na.rm=TRUE, FUN=sd) > 0) )
    x <- x[w,]
    if(!is.null(toporder)) toporder <- toporder[w]
  }
  y <- x
  if(z) y <- t(scale(t(x)))
  if(!is.null(toporder)){
    if(toporder.meth=="before"){
      ag <- aggregate(y, by=list(toporder), na.rm=TRUE, FUN=median)
      row.names(ag) <- ag[,1]
      ag <- ag[,-1]
      if(nrow(ag)>2){
        try( ag <- sortRows(ag, z=FALSE, na.rm=FALSE, method=method),
             silent=TRUE)
      }
      ll <- split(as.data.frame(y), toporder)
      ll <- lapply(ll, FUN=function(x){
          tryCatch(sortRows(x,method=method), error=function(e) return(x))
      })
      y <- unlist(lapply(ll[row.names(ag)], FUN=row.names))
      return(x[y,])
    }else{
      o1 <- get_order(seriate(dist(y), method=method))
      oa <- aggregate(o1,by=list(top=as.factor(toporder)),FUN=median)
      oa <- oa[order(oa[,2]),1]
      toporder <- factor(as.character(toporder), levels=as.character(oa))
      return(x[order(as.numeric(toporder),o1),])
    }
  }
  ss <- seriate(dist(y), method=method)
  x[get_order(ss),]
}


.chooseAssay <- function(se, assayName=NULL){
  if(is.null(assayName) && !is.null(assayNames(se))){
    assayName <- intersect(assayNames(se), c("log2FC", "logFC", "corrected", "imputed", "logcpm", "lognorm"))
    if(length(assayName)>0){
      assayName <- assayName[1]
      message("Using assay ", assayName)
    }else{
      assayName <- NULL
    }
  }
  if(!is.null(assayName) && !any(assayName %in% assayNames(se))) stop("Assay '", assayName, "' not found!")
  if(is.null(assayName)){
    if(length(assays(se))>1) message("Assay unspecified, and multiple assays present - will use the first one.")
    return(assay(se))
  }
  assays(se)[[intersect(assayName,assayNames(se))[1]]]
}

.getHMcols <- function(cols=NULL, n=29){
  if(is.null(cols)) cols <- .getDef("hmcols")
  if(is.function(cols)) return(cols)
  if(length(cols) %in% 2:3)  cols <- colorRampPalette(cols)(n)
  cols
}

#' getBreaks
#'
#' Produces symmetrical breaks for a color scale, with the scale steps
#' increasing for large values, which is useful to avoid outliers influencing
#' too much the color scale.
#'
#' @param x A matrix of log2FC (or any numerical values centered around 0)
#' @param n The desired number of breaks.
#' @param split.prop The proportion of the data points to plot on a linear
#' scale; the remaining will be plotted on a scale with regular frequency per
#' step (quantile).
#'
#' @return A vector of breaks of length = `n`
#' @export
#'
#' @examples
#' dat <- rnorm(100,sd = 10)
#' getBreaks(dat, 10)
getBreaks <- function(x, n, split.prop=0.96){
    x <- abs(x)
    n2 <- floor(n/2)+1
    q <- as.numeric(quantile(x,split.prop))
    xr <- seq(from=0, to=q, length.out=floor(split.prop*n2))
    n2 <- n2-length(xr)
    q <- quantile(as.numeric(x)[which(x>q)],(1:n2)/n2)
    xr <- c(xr,as.numeric(q))
    c(-rev(xr[-1]), xr)
}

.getDef <- function(x){
  a <- c( "Batch", "batch", "Condition","condition", "Group","group",
          "Genotype", "genotype", "Dataset")
  switch(x,
         assay=getOption("SEtools_def_assayName",
                         default=c("logFC", "logcpm", "lognorm")),
         anno_colors=getOption("SEtools_def_anno_colors", default=list()),
         hmcols=getOption("SEtools_def_hmcols",
                          default=c("blue", "black", "yellow")),
         anno_columns=getOption("SEtools_def_anno_columns", default=a),
         anno_rows=getOption("SEtools_def_anno_rows", default=c()),
         gaps_at=getOption("SEtools_def_gaps_at", default="Dataset"),
         breaks=getOption("SEtools_def_breaks", default=NULL)
        )
}


#' resetAllSEtoolsOptions
#'
#' Resents all global options relative to SEtools.
#'
#' @return None
#'
#' @examples
#' resetAllSEtoolsOptions()
#'
#' @export
resetAllSEtoolsOptions <- function(){
  for(o in grep("^SEtools_",names(options()), value=TRUE)){
    eval(parse(text=paste0('options("',o,'"=NULL)')))
  }
}


#' log2FC
#'
#' Generates log2(foldchange) matrix/assay, eventually on a per-batch fashion.
#'
#' @param x A numeric matrix, or a `SummarizedExperiment` object
#' @param fromAssay The assay to use if `x` is a `SummarizedExperiment`
#' @param controls A vector of which samples should be used as controls for
#' foldchange calculations.
#' @param by An optional vector indicating groups/batches by which the controls
#' will be averaged to calculate per-group foldchanges.
#' @param isLog Logical; whether the data is log-transformed. If NULL, will
#' attempt to figure it out from the data and/or assay name
#'
#' @return An object of same class as `x`; if a `SummarizedExperiment`, will
#' have the additional assay `log2FC`.
#'
#' @examples
#' log2FC( matrix(rnorm(40), ncol=4), controls=1:2 )
#'
#' @import SummarizedExperiment
#' @export
log2FC <- function(x, fromAssay=NULL, controls, by=NULL, isLog=NULL){
    if(is.null(colnames(x))) colnames(x) <- paste0("S",seq_len(ncol(x)))
    if(is(x, "SummarizedExperiment")){
        if(is.null(fromAssay))
            stop("If `x` is a SummarizedExperiment, specify the assay to use ",
                "using `fromAssay`")
        if(!(fromAssay %in% assayNames(x)))
            stop("`fromAssay` '", fromAssay, "' not found.")
        if(!is.null(by) && length(by)==1 && by %in% colnames(colData(x)))
            by <- colData(x)[[by]]
        a <- assays(x)[[fromAssay]]
    }else{
        a <- x
    }
    if(is.null(isLog)){
        if(!is.null(fromAssay) && grep("^log",fromAssay, ignore.case=TRUE)){
            isLog <- TRUE
        }else{
            isLog <- any(a<0)
        }
    }
    if(!isLog) a <- log1p(a)
    if(is.logical(controls)) controls <- which(controls)
    if(!all(controls %in% seq_len(ncol(a))))
        stop("Some control indexes are out of range.")
    if(is.null(by)) by <- rep(1,ncol(a))
    i <- split(1:ncol(a),by)
    lfc <- do.call(cbind, lapply(i, FUN=function(x){
        c2 <- intersect(x,controls)
        if(length(c2)==0) stop("Some groups of `by` have no controls.")
        a[,x,drop=FALSE]-rowMeans(a[,c2,drop=FALSE])
    }))
    lfc <- lfc[,colnames(x)]
    if(is(x, "SummarizedExperiment")){
        assays(x)$log2FC <- lfc
        return(x)
    }
    lfc
}
