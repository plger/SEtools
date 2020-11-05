#' svacor
#'
#' A wrapper around SVA-based correction, providing a corrected assay.
#'
#' @param SE An object of class `SummarizedExperiment`. Alternatively, a matrix can be
#' used, but many options will not be supported.
#' @param form The formula of the differential expression model
#' @param form0 An optional formula for the null model
#' @param mm If `form=NULL`, the model.matrix
#' @param mm0 An optional null model.matrix.
#' @param regressOutNull Logical; whether to regress out the variables of `form0` (default TRUE)
#' @param trans Either "svaseq" (default, expects a count assay), "none", or "log".
#' @param assayName Which assay to use (if `SE` is a `SummarizedExperiment`); if missing,
#' will use the first one.
#' @param n.sv The number of surrogate variables (if omitted, `sva` will attempt to
#' estimate it).
#' @param ... Any other argument passed to the `sva` command.
#'
#' @return Returns a `SummarizedExperiment` (with a `corrected` assay and the surrogate
#' variables in `colData`) if `SE` is an object of that class `SummarizedExperiment`;
#' otherwise, a list with the slots:
#' * `sv`: a table of the surrogate variables
#' * `cor`: the corrected data (for plotting)
#' * `mm`: the model.matrix containing, in addition to the specified experimental
#' variables, all detected surrogate variables.
#'
#' @importFrom sva sva svaseq
#' @import SummarizedExperiment
#' @importFrom stats model.matrix
#' @export
svacor <- function(SE, form=NULL, form0=~1, mm=NULL, mm0=NULL, regressOutNull=TRUE,
                   trans=c("svaseq","none", "log"), assayName=NULL, n.sv=NULL, ...){
  trans <- match.arg(trans)
  if( (is.null(form) && is.null(mm)) ||
      (!is.null(form) && !is.null(mm)) ) stop("Only one of `form` or `mm` should be given.")
  if(is.null(mm)){
    if(is(SE,"SummarizedExperiment")){
      CD <- as.data.frame(colData(SE))
      mm <- model.matrix(form, data=CD)
    }else{
      stop("If `form` is used, `SE` should be a SummarizedExperiment.")
    }
  }
  if(is.null(mm0)){
    if(is.null(form0) || !is.null(mm)){
      mm0 <- mm[,1,drop=F]
    }else{
      mm0 <- model.matrix(form0, data=CD)
    }
  }
  if(is(SE, "SummarizedExperiment")){
    if(is.null(assayName)){
      en <- as.matrix(assay(SE))
    }else{
      en <- as.matrix(assays(SE)[[assayName]])
    }
  }else{
    en <- as.matrix(SE)
  }
  if(trans=="log") en <- log2(en+1)

  if(is.null(n.sv) || n.sv>0){
    if(trans=="svaseq"){
      if(any(en<0)) stop("Trying to run 'svaseq', which expects counts, but the data ",
                         "contains negative values.")
      sv <- svaseq(en, mm, mm0, n.sv=n.sv, ...)
    }else{
      sv <- sva(en, mm, mm0, n.sv=n.sv, ...)
    }
    n.sv <- sv$n.sv
    sv <- sv$sv
  }
  if(n.sv==0){
    if(!regressOutNull | ncol(mm0)==1){
      message("Nothing to do!")
      return(NULL)
    }
    X <- as.matrix(mm)
    mm2 <- mm
  }else{
    colnames(sv) <- paste0("SV",1:ncol(sv))
    X <- cbind(mm, sv)
    mm2 <- cbind(mm[,1,drop=F],sv,mm[,-1,drop=F])
  }
  H <- solve(t(X)%*%X)%*%t(X)
  b <- (H%*%t(en))
  if(regressOutNull){
    cn <- setdiff(colnames(X),setdiff(colnames(mm), colnames(mm0)))
  }else{
    cn <- setdiff(colnames(X),colnames(mm))
  }
  cn <- setdiff(cn, "(Intercept)")
  encor <- en - t(as.matrix(X[,cn]) %*% b[cn,])
  if(is(SE, "SummarizedExperiment")){
    SE <- SE[row.names(encor),]
    if(n.sv>0) colData(SE) <- cbind(colData(SE), sv)
    assays(SE)$corrected <- encor
    return(SE)
  }
  mm2 <- cbind(mm[,1,drop=F],sv,mm[,-1,drop=F])
  return(list(sv=sv, cor=encor, mm=mm2))
}
