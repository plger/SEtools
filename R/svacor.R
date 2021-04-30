#' svacor
#'
#' A wrapper around SVA-based correction, providing a corrected assay. If this is RNAseq
#' data or similar, use a count assay assay with `useVST=TRUE`; otherwise (e.g.
#' proteomics) a log-normalized assay is recommended.
#'
#' @param SE An object of class `SummarizedExperiment`.
#' @param form The formula of the differential expression model
#' @param form0 An optional formula for the null model
#' @param assayName The name (or index) of the assay to use.
#' @param regressOutNull Logical; whether to regress out the variables of `form0`.
#' @param useVST Logical; whether to use DESeq2's variance-stabilizing transformation;
#' (for count data!)
#' @param n.sv The number of surrogate variables (if omitted, \code{\link{sva}} will
#' attempt to estimate it)
#' @param ... Any other argument passed to the \code{\link{sva}} command.
#'
#' @return Returns the `SummarizedExperiment` with a `corrrected` assay and the surrogate
#' variables in `colData`.
#'
#' @importFrom DESeq2 DESeqDataSetFromMatrix estimateSizeFactors vst varianceStabilizingTransformation
#' @importFrom sva sva
#' @importFrom stats model.matrix
#' @import SummarizedExperiment
#' @export
svacor <- function(SE, form, form0=~1, assayName=NULL, regressOutNull=TRUE, useVST=TRUE,
                   n.sv=NULL, ...){
  if(!is(SE,"SummarizedExperiment")) stop("`SE` should be a SummarizedExperiment.")
  CD <- as.data.frame(colData(SE))
  mm <- model.matrix(form, data=CD)
  mm0 <- model.matrix(form0, data=CD)

  if(is.null(assayName)){
    if(useVST && any(assayNames(SE)=="counts")){
      assayName <- "counts"
    }else{
      message("assayName not specified, using the first available.")
      assayName <- 1
    }
  }
  en <- as.matrix(assays(SE)[[assayName]])

  if(useVST){
    message("Using variance-stabilizing transformation")
    en <- tryCatch({
      dds <- DESeqDataSetFromMatrix(round(en), CD, form)
      dds <- estimateSizeFactors(dds)
      as.matrix(assay(vst(dds, blind=FALSE)))
    }, error=function(e){
      varianceStabilizingTransformation(round(en))
    })
  }
  if(is.null(n.sv) || n.sv>0){
    sv <- sva(en, mm, mm0, n.sv=n.sv, ...)
    n.sv <- sv$n.sv
    sv <- sv$sv
  }
  if(n.sv==0){
    if(!regressOutNull | ncol(mm0)==1){
      message("Nothing to do!")
      return(SE)
    }
    X <- as.matrix(mm)
    mm2 <- mm
  }else{
    colnames(sv) <- paste0("SV",1:ncol(sv))
    X <- cbind(mm, sv)
    mm2 <- cbind(mm[,1,drop=FALSE],sv,mm[,-1,drop=FALSE])
  }
  H <- solve(t(X) %*% X) %*% t(X)
  b <- (H %*% t(en))
  if(regressOutNull){
    cn <- setdiff(colnames(X),setdiff(colnames(mm), colnames(mm0)))
  }else{
    cn <- setdiff(colnames(X),colnames(mm))
  }
  cn <- setdiff(cn, "(Intercept)")
  encor <- en - t(as.matrix(X[,cn]) %*% b[cn,])
  SE <- SE[row.names(encor),]
  if(length(pSVs <- grep("^SV[0-9]+$", colnames(colData(SE)), value=TRUE))>0){
    warning("Found and removed previous SV columns in colData.")
    colData(SE) <- colData(SE)[,setdiff(colnames(colData(SE)), pSVs),drop=FALSE]
  }
  if(n.sv>0) colData(SE) <- cbind(colData(SE), sv)
  assays(SE)$corrected <- encor
  SE
}
