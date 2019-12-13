#' mergeSEs
#'
#' Merges a list of
#' \code{\link[SummarizedExperiment]{SummarizedExperiment-class}}, either by
#' row.names or through specified rowData fields. In cases of many-to-many 
#' (or one-to-many) mappings, `aggFun` determines whether the records are 
#' aggregated by linking ID (if an aggregation method is given) or all 
#' combinations are returned (if `aggFun=NULL` - default).
#'
#' @param ll A (named) list of
#' \code{\link[SummarizedExperiment]{SummarizedExperiment-class}}
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
#' @param mergeBy The `rowData` column to merge with. If NULL, row.names are
#' used.
#' @param aggFun The aggregation function to use when multiple rows have the
#' same `mergeBy` value. If merging multiple assays, a different function per
#' assay can be passed as a named list (see \code{\link[SEtools]{aggSE}}). If
#' NULL (default), entries will be reused to have each combination.
#' @param addDatasetPrefix Logical; whether the name of the dataset should be
#' appended to the sample names (default TRUE).
#' @param defValues An optional named list of default `colData` values when some
#'  columns are missing from some SEs.
#'
#' @return An object of class
#' \code{\link[SummarizedExperiment]{SummarizedExperiment-class}}
#'
#' @examples
#' data("SE", package="SEtools")
#' mergeSEs( list( se1=SE[,1:10], se2=SE[,11:20] ) )
#'
#' @import SummarizedExperiment BiocParallel
#' @importFrom data.table data.table rbindlist
#' @export
mergeSEs <- function(ll, use.assays=NULL, do.scale=TRUE, commonOnly=TRUE,
                     colColumns=NULL, mergeBy=NULL, aggFun=NULL,
                     addDatasetPrefix=TRUE, defValues=list(),
                     BPPARAM=SerialParam()){
  ll <- .forceAssayNames(ll)
  if(!is.null(mergeBy)){
      rdf <- table(unlist(lapply(ll, FUN=function(x) colnames(rowData(x)))))
      if(!all(mergeBy %in% names(rdf)) || !all(rdf[mergeBy]==length(ll)))
          stop("The `mergeBy` field(s) are missing from the rowData of some
               object(s).")
      if(is.null(aggFun)){
          se <- .DFmerge_mult(ll, fields=mergeBy, use.assays=use.assays,
                              all=!commonOnly)
      }else{
          ll <- lapply(ll, an=.commonAssays(ll, use.assays),
                       FUN=function(x,an){
                           assays(x) <- assays(x)[an]
                           .has_nan
                           x
                       })
          message("Aggregating the objects by ", paste(mergeBy,collapse=", "))
          ll <- bplapply(ll, by=mergeBy, assayFun=aggFun, BPPARAM=BPPARAM,
                         FUN=aggSE)
          message("Merging...")
          mergeBy <- NULL
      }
  }
  if(is.null(mergeBy)){
      tt <- table(unlist(lapply(ll,row.names)))
      if(max(tt)==1) stop("No matching row.names!")
      if(commonOnly){
          g <- names(tt)[which(tt==length(ll))]
      }else{
          g <- names(tt)
      }
      dat <- .prepAssays(ll, use.assays=use.assays, do.scale=do.scale, rn=g)
      cd2 <- .mergeColData(ll, colColumns=colColumns, defValues=defValues)
      se <- SummarizedExperiment(dat, rowData=.mergeRowData(ll, g), colData=cd2)
  }

  metadata(se) <- .mergeMetadata(ll)

  nn <- do.call(c, lapply(ll, colnames))
  if(!is.null(names(ll)) &&
     ( addDatasetPrefix || any(table(nn)>1) ) ){
    nn <- paste(rep(names(ll), sapply(ll,ncol)), nn, sep=".")
  }
  colnames(se) <- nn
  se
}

.forceAssayNames <- function(ll){
    if(any(sapply(lapply(ll,assayNames),is.null)))
        warning("Some objects have no assay names - arbitrary names are given.")
    lapply(ll, FUN=function(x){
        if(is.null(assayNames(x)))
            assayNames(x) <- paste0("assay", seq_len(length(assays(x))))
        x
    })
}

.commonAssays <- function(ll, use.assays=NULL){
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
       warning( "Some of the given assays were not available and will
                 be discarded: ",
                 paste(setdiff(old, use.assays), collapse=", ") )
    }
  }
  if(length(use.assays)==0){
    if(is.null(names(ll))) names(ll) <- paste("object", seq_len(length(ll)))
      aa <- paste0(names(ll), ":\n",
                   sapply(ll, FUN=function(x) paste(assayNames(x), collapse=", ")),
                   collapse="\n" )
      stop("No assay to merge. Available assays:\n", aa, "\n",
           "To merge the first assay of each object, use `use.assays=1`.")
  }
  use.assays
}

.prepAssays <- function(ll, use.assays=NULL, do.scale=FALSE, rn=NULL){
  ll <- .forceAssayNames(ll)
  use.assays <- .commonAssays(ll, use.assays=use.assays)
  ll <- lapply(ll, FUN=function(x){
      assays(x) <- assays(x)[use.assays]
      x
  })
  if(is.null(rn)) rn <- unique(unlist(sapply(ll, FUN=row.names)))
  if(length(do.scale)==1) do.scale <- rep(do.scale, length(use.assays))
  if(length(do.scale)!=length(use.assays))
    stop( "`do.scale` should have a length either of 1 or equal to the number ",
          "of assays used.")

  a <- lapply(seq_len(length(use.assays)), FUN=function(a){
    x <- lapply(ll, FUN=function(x){
      x <- assays(x)[[use.assays[[a]]]]
      if(all(rn %in% row.names(x))) return(x[rn,])
      as.matrix(as.data.frame(x[rn,,drop=FALSE]))
    })
    if(do.scale[a]){
      if(any(sapply(x, FUN=function(x) any(is.infinite(x) | is.na(x))))){
        stop("Cannot scale the data in the presence of missing or infinite
             values.")
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

.mergeMetadata <- function(ll){
  if(is(ll[[1]], "SummarizedExperiment")) ll <- lapply(ll, metadata)
  out <- lapply(ll, FUN=function(x) x[setdiff(names(x),"anno_colors")])
  x <- lapply(ll, FUN=function(x) x[["anno_colors"]])
  x <- x[!sapply(x,is.null)]
  out$anno_colors <- .mergelists(x)
  out
}

.mergeColData <- function(ll, colColumns=NULL, defValues=list()){
    cd <- lapply(ll, cc=colColumns, FUN=function(x,cc){
        x <- as.data.frame(colData(x))
        if(!is.null(cc)) x <- x[,intersect(cc,colnames(x)),drop=FALSE]
        x
    })
    cd2 <- as.data.frame(rbindlist(cd, fill=TRUE, idcol="Dataset"))
    for(i in names(defValues))
        cd2[[i]][which(is.na(cd2[[i]]))] <- defValues[[i]]
    cd2
}

.mergeRowData <- function(ll, g){
  rd <- lapply(ll, FUN=function(x){
    if(is(x,"SummarizedExperiment")) x <- rowData(x)
    if(!is.null(g)) x <- as.data.frame(x)[g,,drop=FALSE]
    x
  })
  if(is.null(names(rd))) names(rd) <- paste0("D",seq_len(length(rd)))
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
    w <- which(sapply(seq_along(x), FUN=function(i){
      any(sapply(seq_len(i-1), FUN=function(j,i){
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

.DFmerge_mult <- function(SEs, fields=NULL, use.assays=NULL, do.scale=FALSE,
                          colColumns=NULL, defValues=list(), all=FALSE){
    use.assays <- .commonAssays(SEs, use.assays=use.assays)
    ll <- lapply(SEs, rowData)
    if(is.null(fields)){
        fields <- table(unlist(lapply(ll,colnames)))
        fields <- names(fields)[which(fields==length(ll))]
        if(length(fields)==0) stop("The objects share no rowData columns.")
        message("Merging by ", paste(fields, collapse=", "))
    }
    if(!all(vapply(ll, FUN.VALUE=logical(1), FUN=function(x)
        all(fields %in% colnames(x)) )))
            stop("Some fields are not in all objects.")
    ll <- lapply(names(ll), FUN=function(n){
        x <- ll[[n]]
        x[[paste0(n,".id")]] <- row.names(x)
        x
    })
    m <- suppressWarnings({
        Reduce(function(x, y) merge(x, y, by=fields, all=all), ll)
    })
    m <- m[,!duplicated(colnames(m))]
    fields <- c(fields, grep("\\.id$", colnames(m), value=TRUE))
    row.names(m) <- do.call(paste, m[,fields,drop=FALSE])
    row.names(m) <- sapply( strsplit(row.names(m)," |\\."), FUN=function(x){
        paste(unique(x), collapse=".")
    })

    cd2 <- .mergeColData(SEs, colColumns=colColumns, defValues=defValues )

    a <- lapply(use.assays, FUN=function(an){
        do.call(cbind, lapply(names(SEs), FUN=function(n){
            x <- .prepAssays(SEs[n], use.assays=an, do.scale=do.scale)[[1]]
            x <- as.matrix(as.data.frame(x)[m[[paste0(n,".id")]],])
            row.names(x) <- row.names(m)
            x
        }))
    })
    names(a) <- use.assays
    se <- SummarizedExperiment( a, colData=cd2, rowData=m )
    se
}

