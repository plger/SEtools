#' sortRows
#'
#' @param x A numeric matrix or data.frame.
#' @param z Whether to scale rows for the purpose of calculating order (default FALSE).
#' @param toporder Optional verctor of categories (length=nrow(x)) on which to supra-order
#'  when sorting rows.
#' @param na.rm Wheter to remove missing values and invariant rows (default FALSE).
#' @param method Seriation method; 'MDS_angle' (default) or 'R2E' recommended.
#' @param toporder.meth Whether to perform higher-order sorting 'before' (default) or
#' 'after' the lower-order sorting.
#'
#' @return A reordered matrix or data.frame.
#'
#' @importFrom seriation seriate get_order
#' @export
sortRows <- function(x, z=F, toporder=NULL, na.rm=F, method="MDS_angle", toporder.meth="before"){
  toporder.meth <- match.arg(toporder.meth, c("before","after"))
  if(na.rm){
    w <- which( apply(x, 1, FUN = function(y){ !any(is.na(y)) }) |
                  !(apply(x, 1, na.rm=T, FUN=sd) > 0) )
    x <- x[w,]
    if(!is.null(toporder)) toporder <- toporder[w]
  }
  y <- x
  if(z) y <- t(scale(t(x)))
  if(!is.null(toporder)){
    if(toporder.meth=="before"){
      ag <- aggregate(y, by=list(toporder), na.rm=T, FUN=median)
      row.names(ag) <- ag[,1]
      ag <- sortRows(ag[,-1], z=F, na.rm=F, method=method)
      ll <- split(as.data.frame(y), toporder)
      ll <- lapply(ll, z=F, na.rm=F, method=method, FUN=sortRows)
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


