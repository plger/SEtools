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
      ag <- sortRows(ag[,-1], z=FALSE, na.rm=FALSE, method=method)
      ll <- split(as.data.frame(y), toporder)
      ll <- lapply(ll, z=FALSE, na.rm=FALSE, method=method, FUN=sortRows)
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

.getBreaks <- function(x, n){
  if(!is.list(x)) x <- list(x)
  xr <- range(sapply(x, na.rm=TRUE, FUN=range), na.rm=TRUE)
  if(ceiling(max(abs(xr)))<=2){
    xr <- ceiling(max(abs(xr*10)))/10
  }else{
    xr <- ceiling(max(abs(xr)))
  }
  if(xr>=4){
    breaks <- c( -xr, -3.5, -3,
                 seq(from=-2.5,to=2.5,length.out=n-7),
                 3, 3.5, xr)
  }else{
    if(xr>=3){
      breaks <- c(-xr, -2.5,
                  seq(from=-2,to=2,length.out=n-5),
                  2.5, xr)
    }else{
      breaks <- seq(from=-xr,to=xr,length.out=n-1)
    }
  }
  breaks
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
#' @return NULL
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

