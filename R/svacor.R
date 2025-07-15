#' svacor
#'
#' A wrapper around SVA-based correction, providing a corrected assay. If this
#' is RNAseq data or similar, use a count assay assay with `method` either 'vst'
#' or 'svaseq'; otherwise (e.g. proteomics) a log-normalized assay is
#' recommended with `method="sva"`. Note that the corrected assay, while useful
#' for visualization, should be interpreted with care, as they omit major
#' variation!
#'
#' @param SE An object of class `SummarizedExperiment`.
#' @param form The formula of the differential expression model
#' @param form0 An optional formula for the null model
#' @param assayName The name (or index) of the assay to use.
#' @param regressOutNull Logical; whether to regress out the variables of `form0`.
#' @param useVST Deprecated; use the `method` argument instead.
#' @param method Either 'vst' (uses DESeq2 variance-stabilization before
#'   running SVA), 'svaseq' (uses sva::svaseq on normalized counts), or 'sva'
#'   (uses standard sva, not appropriate for count data).
#' @param n.sv The number of surrogate variables (if omitted, \code{\link{sva}}
#' will attempt to estimate it). Note that automatic determination of the number
#' of SVs will often lead to fairly large number of SVs, use
#' `numSVmethod="leek"` for a more conservative estimate.
#' @param ... Any other argument passed to the \code{\link{sva}} command.
#'
#' @return Returns the `SummarizedExperiment` with a `corrrected` assay and the
#' surrogate variables in `colData`.
#'
#' @importFrom DESeq2 DESeqDataSetFromMatrix estimateSizeFactors vst
#' @importFrom DESeq2 varianceStabilizingTransformation
#' @importFrom sva sva svaseq
#' @importFrom stats model.matrix
#' @import SummarizedExperiment
#' @export
#' @examples
#' data("SE", package="SEtools")
#' SE <- svacor(SE, ~Condition)
svacor <- function(SE, form, form0=~1, assayName=NULL, regressOutNull=TRUE,
                   method=c("vst","svaseq","sva"), useVST=NULL, n.sv=NULL, ...){
  if(!is(SE,"SummarizedExperiment")) stop("`SE` should be a SummarizedExperiment.")
  if(!is.null(useVST)){
      stop("The `useVST` argument is deprecated, please use the `method` ",
           "argument instead.")
  }
  method <- match.arg(method)
  CD <- as.data.frame(colData(SE))
  mm <- model.matrix(form, data=CD)
  mm0 <- model.matrix(form0, data=CD)

  if(is.null(assayName)){
    if(method %in% c("vst","svaseq") && any(assayNames(SE)=="counts")){
      assayName <- "counts"
    }else{
      message("assayName not specified, using the first available.")
      assayName <- 1
    }
  }
  en <- as.matrix(assays(SE)[[assayName]])

  if(method=="vst"){
    message("Using variance-stabilizing transformation")
    en <- tryCatch({
      dds <- DESeqDataSetFromMatrix(round(en), CD, form)
      dds <- estimateSizeFactors(dds)
      as.matrix(assay(vst(dds, blind=FALSE)))
    }, error=function(e){
      varianceStabilizingTransformation(round(en))
    })
  }else if(method=="svaseq"){
    en <- .normcounts(en)
  }
  if(is.null(n.sv) || n.sv>0){
    if(method=="svaseq"){
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

#' @importFrom edgeR calcNormFactors cpm DGEList
.normcounts <- function(x){
    if(is(x,"SummarizedExperiment")) x <- assay(x)
    d <- calcNormFactors(DGEList(x))
    return( edgeR::cpm(d, normalized.lib.sizes=TRUE) *
                median(d$samples$lib.size)/1000000 )
}
