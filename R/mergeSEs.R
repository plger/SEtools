#' mergeSEs
#'
#' Merges a list of `SummarizedExperiment`.
#'
#' @param ll A (named) list of SummarizedExperiments
#' @param use.assays Names (or indexes) of the assays to use. By default, all
#' common assays are used.
#' @param do.scale A logical vector indicating (globally or for each assay)
#' whether to perform row unit-variance scaling on each dataset before merging
#' (default TRUE).
#' @param commonOnly Logical; whether to restrict to rows present in all
#' datasets (default TRUE).
#' @param colColumns A character vector specifying `colData` columns to include
#' (if available in at least one of the datasets). If NULL, everything is kept.
#' @param defValues A list specifying the default values for `colColumns` when
#' these are absent.
#' @param addDatasetPrefix Logical; whether the name of the dataset should be
#' appended to the sample names (default TRUE).
#'
#' @return An object of class `SummarizedExperiment`
#'
#' @examples
#' data("SE", package="SEtools")
#' mergeSEs( list( se1=SE[,1:10], se2=SE[,11:20] ) )
#'
#' @import SummarizedExperiment
#' @importFrom data.table data.table rbindlist
#' @export
mergeSEs <- function(ll, use.assays=NULL, do.scale=TRUE, commonOnly=TRUE,
                     colColumns=NULL, addDatasetPrefix=TRUE, defValues=list()){
  if(!commonOnly && any(do.scale)) stop("For z-scores, `commonOnly` must be TRUE.")

  tt <- table(unlist(lapply(ll,row.names)))
  if(max(tt)==1) stop("No matching row.names!")
  if(commonOnly){
    g <- names(tt)[which(tt==length(ll))]
  }else{
    g <- names(tt)
  }

  dat <- .prepAssays(ll, use.assays=use.assays, do.scale=do.scale, rn=g)

  # rowData
  rd <- .mergeRowData(ll, g)

  # colData
  cd <- lapply(ll, cc=colColumns, FUN=function(x,cc){
    x <- as.data.frame(colData(x))
    if(!is.null(cc)) x <- x[,intersect(cc,colnames(x)),drop=FALSE]
    x
  })
  cd2 <- as.data.frame(data.table::rbindlist(cd, fill=T, idcol="Dataset"))
  for(i in names(defValues)) cd2[[i]][which(is.na(cd2[[i]]))] <- defValues[[i]]

  se <- SummarizedExperiment( dat, colData=cd2, rowData=rd )
  nn <- do.call(c, lapply(ll, colnames))
  if(!is.null(names(ll)) &&
     ( addDatasetPrefix || any(table(nn)>1) ) ){
    nn <- paste(rep(names(ll), sapply(ll,ncol)), nn, sep=".")
  }
  colnames(se) <- nn
  se
}



.prepAssays <- function(ll, use.assays=NULL, do.scale=FALSE, rn=NULL){
  if(is.null(rn)) rn <- unique(unlist(sapply(ll, FUN=row.names)))
  an <- table(unlist(lapply(ll, FUN=assayNames)))
  if(is.null(use.assays)){
    use.assays <- names(an)[which(an==length(ll))]
  }
  if(is.numeric(use.assays)){
    if(max(use.assays) > min(sapply(ll, FUN=function(x) length(assays(x))))){
      stop("Some assays requested are not available in some of the objects.")
    }
  }else{
    an <- names(an)[which(an==length(ll))]
    old <- use.assays
    use.assays <- intersect(use.assays, an)
    if(length(old)!=length(use.assays)){
      warning( "Some of the given assays were not available and will be discarded: ",
               paste(setdiff(old, use.assays), collapse=", ") )
    }
  }
  if(length(use.assays)==0){
    if(is.null(names(ll))) names(ll) <- paste("object", 1:length(ll))
    aa <- paste0(names(ll),":\n", sapply(ll, FUN=function(x) paste(assayNames(x), collapse=", ")), collapse="\n")
    stop("No assay to merge. Available assays:\n", aa, "\n",
         "To merge the first assay of each object, use `use.assays=1`.")
  }
  if(length(do.scale)==1) do.scale <- rep(do.scale, length(use.assays))
  if(length(do.scale)!=length(use.assays))
    stop("`do.scale` should have a length either of 1 or equal to the number of assays used.")

  a <- lapply(1:length(use.assays), FUN=function(a){
    x <- lapply(ll, FUN=function(x){
      x <- assays(x)[[use.assays[[a]]]]
      if(all(rn %in% row.names(x))) return(x[rn,])
      as.matrix(as.data.frame(x[rn,,drop=FALSE]))
    })
    if(do.scale[a]){
      if(any(sapply(x, FUN=function(x) any(is.infinite(x) | is.na(x))))){
        stop("Cannot scale the data in the presence of missing or infinite values.")
      }
      x <- tryCatch(  lapply(x, FUN=function(x) t(scale(t(x))) ),
                      error=function(e){
                        lapply(x, FUN=function(x){
                          apply(x,1,FUN=function(x){
                            if(!(sd(x)>0)) return(rep(0,length(x)))
                            as.numeric(scale(as.numeric(x)))
                          })
                        })
                      })
    }
    do.call(cbind, x)
  })
  if(!is.numeric(use.assays)) names(a) <- use.assays
  a
}

.mergeRowData <- function(ll, g){
  rd <- lapply(ll, FUN=function(x){
    as.data.frame(rowData(x))[g,,drop=FALSE]
  })
  if(is.null(names(rd))) names(rd) <- paste0("D",1:length(rd))
  cn <- unique(unlist(lapply(rd, FUN=colnames)))
  rd <- lapply(cn, FUN=function(n){ # for each unique rowData column name
    x <- lapply(rd, FUN=function(y){
      if(n %in% colnames(y)) return(y[,n])
      return(NULL)
    })
    x <- x[which(!sapply(x,is.null))]
    if(length(x)==0) return(NULL)
    if(length(x)==1){
      x <- data.frame(x[[1]])
      colnames(x) <- n
      return(x)
    }

    # the column occurs in more than one dataset
    # we first remove identical columns
    w <- which(sapply(2:length(x), FUN=function(i){
      any(sapply(1:(i-1), FUN=function(j,i){
        all(x[[i]]==x[[j]])
      }))
    }))

    coln <- paste(names(x), n, sep=".")
    if(length(w)>0){
      x <- x[-w]
      if(length(x)==1){
        x <- data.frame(x[[1]])
        colnames(x) <- n
        return(x)
      }
      coln <- coln[-w]
    }
    x <- do.call(cbind, x)
    colnames(x) <- coln
    x
  })
  rd <- rd[!sapply(rd,is.null)]
  if(length(rd)==0) return(NULL)
  rd <- do.call(cbind, rd)
  row.names(rd) <- g
  rd
}
