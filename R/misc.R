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
  if(is.factor(toporder)) toporder <- droplevels(toporder)
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


.chooseAssay <- function(se, assayName=NULL, returnName=FALSE){
  a <- .getDef("assay")
  if(is.null(assayName) && !is.null(assayNames(se))){
    assayName <- intersect(assayNames(se), a)
    if(length(assayName)>0){
      assayName <- assayName[1]
      message("Using assay ", assayName)
    }else{
      assayName <- NULL
    }
  }
  if(!is.null(assayName) && !is.numeric(assayName) &&
     !any(assayName %in% assayNames(se)))
      stop("Assay '", assayName, "' not found!")
  if(is.null(assayName)){
    if(length(assays(se))>1) message("Assay unspecified, and multiple assays",
                                        " present - will use the first one.")
    assayName <- 1
  }else{
    assayName <- intersect(assayName,assayNames(se))[1]
  }
  if(returnName) return(assayName)
  assays(se)[[assayName]]
}

.getHMcols <- function(cols=NULL, n=100){
  if(is.null(cols)) cols <- .getDef("hmcols")
  if(is.function(cols)) return(cols)
  if(length(cols) %in% 2:3)  return(colorRampPalette(cols)(n))
  cols
}

.getBaseHMcols <- function(se, cols){
    if(!is.null(cols)) return(cols)
    if(!is.null(se) && !is.null(cols <- metadata(se)$hmcols)) return(cols)
    .getDef("hmcols")
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
#' @param symmetric Logical; whether breaks should be symmetric around 0
#'  (default TRUE)
#'
#' @return A vector of breaks of length = `n`
#' @export
#'
#' @examples
#' dat <- rnorm(100,sd = 10)
#' getBreaks(dat, 10)
getBreaks <- function(x, n, split.prop=0.98, symmetric=TRUE){
    if(is.logical(split.prop)) split.prop <- ifelse(split.prop,0.98,1)
    if(symmetric){
        x <- abs(x)
        n2 <- floor(n/2)+1
    }else{
        n2 <- n
    }
    q <- as.numeric(quantile(x,split.prop,na.rm=TRUE))
    xr <- seq(from=0, to=q, length.out=floor(split.prop*n2))
    n2 <- n2-length(xr)
    if(n2>0){
        q <- quantile(as.numeric(x)[which(x>q)],(1:n2)/n2, na.rm=TRUE)
        xr <- c(xr,as.numeric(q))
    }
    if(symmetric) xr <- c(-rev(xr[-1]), xr)
    if(any(duplicated(xr))){
        ## duplicated breaks, probably because we have to few datapoints;
        ## we fall back onto a linear scale
        xr <- getBreaks(x, n, 1, symmetric=symmetric)
    }
    xr
}

.getDef <- function(x, se){
  a <- c( "Batch", "batch", "Condition","condition", "Group", "group", "Dataset",
          "Genotype", "genotype", "cluster_id", "group_id", "celltype")
  switch(x,
         assay=getOption("SEtools_def_assayName",
                         default=c("logFC", "log2FC", "logcpm", "lognorm")),
         anno_colors=getOption("SEtools_def_anno_colors", default=list()),
         hmcols=getOption("SEtools_def_hmcols",
                          default=c("blue", "black", "yellow")),
         anno_columns=getOption("SEtools_def_anno_columns", default=a),
         anno_rows=getOption("SEtools_def_anno_rows", default=c()),
         gaps_at=getOption("SEtools_def_gaps_at",
                           default=c("Dataset","cluster_id")),
         breaks=getOption("SEtools_def_breaks", default=NULL)
        )
}

.getAnnoCols <- function(se, given=list(), do.assign=FALSE){
    ll <- list( default=.getDef("anno_colors") )
    if(!is.null(metadata(se)$anno_colors)) ll$object <- metadata(se)$anno_colors
    ll$given <- given
    ac <- .mergelists(ll)
    if(do.assign) ac <- .assignAnnoColors(se, ac)
    lapply(ac, unlist)
}

#' @importFrom randomcoloR distinctColorPalette
#' @importFrom SummarizedExperiment colData rowData
.assignAnnoColors <- function(x, anno_colors){
    fn <- function(x){
        if(is.factor(x)) return(levels(x))
        if(is.character(x)) return(unique(x))
        return(NULL)
    }
    if(is(x, "SummarizedExperiment")){
        df <- c( lapply(colData(x), fn), lapply(rowData(x), fn) )
    }else{
        df <- lapply( x, fn)
    }
    for(f in names(df)){
        if(!(f %in% names(anno_colors))) anno_colors[[f]] <- list()
        x <- setdiff(df[[f]], anno_colors[[f]])
        if(length(x)>0) anno_colors[[f]][x] <- distinctColorPalette(length(x))
    }
    anno_colors
}


# non recursive, latest values win
.mergelists <- function(ll){
    names(ll) <- NULL
    names(nn) <- nn <- unique(unlist(lapply(ll,names)))
    lapply(nn, FUN=function(x){
        x <- lapply(ll, function(y) y[[x]])
        x <- do.call(c,x)
        x[!duplicated(names(x))]
    })
}

.has_nan <- function(x){
    if(is(x,"SummarizedExperiment"))
        return(any( sapply(assays(x), .has_nan) ))
    any(is.infinite(x) | is.na(x))
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
#' @param agFun Aggregation function for the baseline (default rowMeans)
#' @param toAssay The name of the assay in which to save the output.
#'
#' @return An object of same class as `x`; if a `SummarizedExperiment`, will
#' have the additional assay named from `toAssay`.
#'
#' @examples
#' log2FC( matrix(rnorm(40), ncol=4), controls=1:2 )
#'
#' @import SummarizedExperiment
#' @export
log2FC <- function(x, fromAssay=NULL, controls, by=NULL, isLog=NULL,
                   agFun=rowMeans, toAssay="log2FC"){
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
	if(!is.matrix(x))
	  stop("`x` should either be a SummarizedExperiment or a numeric matrix.")
        a <- x
    }
    if(is.null(isLog)){
        if(!is.null(fromAssay) && grep("^log",fromAssay, ignore.case=TRUE)){
            isLog <- TRUE
        }else{
            isLog <- any(a<0)
        }
    }
    if(!isLog) a <- log2(a+1)
    if(is.logical(controls)) controls <- which(controls)
    if(!all(controls %in% seq_len(ncol(a))))
        stop("Some control indexes are out of range.")
    if(is.null(by)) by <- rep(1,ncol(a))
    i <- split(1:ncol(a),by)
    lfc <- do.call(cbind, lapply(i, FUN=function(x){
        c2 <- intersect(x,controls)
        if(length(c2)==0) stop("Some groups of `by` have no controls.")
        a[,x,drop=FALSE]-agFun(a[,c2,drop=FALSE],na.rm=TRUE)
    }))
    lfc <- lfc[,colnames(x)]
    if(is(x, "SummarizedExperiment")){
        assays(x)[[toAssay]] <- lfc
        return(x)
    }
    lfc
}

#' flattenPB
#'
#' Flattens a pseudo-bulk SummarizedExperiment as produced by
#' `muscat::aggregateData` so that all cell types are represented in a single
#' assay. Optionally normalizes the data and calculates per-sample logFCs.
#'
#' @param pb a pseudo-bulk SummarizedExperiment as produced by
#' `muscat::aggregateData`, with different celltypes/clusters are assays.
#' @param getLFC Logical; whether to compute `logcpm` and `log2FC` assays.
#'
#' @return A SummarizedExperiment
#' @importFrom edgeR cpm calcNormFactors DGEList
#' @import SummarizedExperiment S4Vectors
#' @export
flattenPB <- function(pb, norm=TRUE, lfc_group="group_id"){
    a <- do.call(cbind, as.list(assays(pb)))
    v.samples <- rep(colnames(pb),length(assays(pb)))
    v.clusters <- rep(assayNames(pb),each=ncol(pb))
    colnames(a) <- paste( v.samples, v.clusters, sep="." )
    cd <- do.call(rbind, lapply(seq_along(assays(pb)),
                                FUN=function(x) colData(pb)) )
    row.names(cd) <- colnames(a)
    cd$cluster_id <- v.clusters
    se <- SummarizedExperiment( list(counts=a), colData=cd, rowData=rowData(pb))
    se$metadata <- pb$metadata
    if(!is.null(metadata(pb)$n_cells)){
        n_cells <- tryCatch({
            mapply( as.character(v.clusters), as.character(v.samples),
                    FUN=function(x,y) metadata(pb)$n_cells[x,y] )
        }, error=function(e){ warning(e); NULL} )
        if(!is.null(n_cells)) se$n_cells <- as.numeric(n_cells)
    }
    if(norm) assays(se)$logcpm <-
        log2(edgeR::cpm(calcNormFactors(DGEList(assay(se))))+1)
    if(is.null(lfc_group) || is.na(lfc_group)) return(se)
    if(is.null(se[[lfc_group]])){
        warning("Could not find '",lfc_group,"', and did not compute log2FC assay.")
        return(se)
    }
    if(!is.factor(se[[lfc_group]])){
        se[[lfc_group]] <- factor(se[[lfc_group]])
        message("Using '", levels(se[[lfc_group]])[1],
                "' as baseline condition")
    }
    log2FC(se, "logcpm", se[[lfc_group]]==levels(se[[lfc_group]])[1],
                 by=se$cluster_id)
}


#' se2xlsx
#'
#' Writes a SummarizedExperiment to an excel/xlsx file. Requires the `openxlsx`
#' package.
#'
#' @param se The `SummarizedExperiment`
#' @param filename 	xlsx file name
#' @param addSheets An optional list of additional tables to save as sheets.
#'
#' @return Saves to file.
#'
#' @examples
#' data("SE", package="SEtools")
#' # not run
#' # se2xls(SE, filename="SE.xlsx")
#'
#' @export
se2xls <- function(se, filename, addSheets=NULL){
    library(openxlsx)
    a <- list( sample_annotation=as.data.frame(colData(se)) )
    if(ncol(rowData(se))>0) a$feature_annotation=as.data.frame(rowData(se))
    a <- c(a, as.list(assays(se)), addSheets)
    write.xlsx(a, file=filename, row.names=TRUE, col.names=TRUE)
}


.prepareAnnoDF <- function(an, anno_colors, fields, whichComplex=NULL,
                           show_legend=TRUE, show_annotation_name=TRUE,
                           dropEmptyLevels=TRUE){
    if(!is.null(whichComplex))
        whichComplex <- match.arg(whichComplex, c("row","column"))
    an <- as.data.frame(an)
    an <- an[,intersect(fields, colnames(an)),drop=FALSE]
    if(ncol(an)==0){
        an <- NULL
    }else{
        for(i in colnames(an)){
            if(is.factor(an[[i]])){
                if(dropEmptyLevels) an[[i]] <- droplevels(an[[i]])
            }
            if(is.logical(an[[i]])){
                an[[i]] <- factor(as.character(an[[i]]),levels=c("FALSE","TRUE"))
                if(!(i %in% names(anno_colors))){
                    anno_colors[[i]] <- c("FALSE"="white", "TRUE"="darkblue")
                }
            }else{
                if(i %in% names(anno_colors)){
                    w <- intersect(names(anno_colors[[i]]),unique(an[[i]]))
                    if(length(w)==0){
                        anno_colors[[i]] <- NULL
                    }else{
                        anno_colors[[i]] <- anno_colors[[i]][w]
                    }
                }
            }
        }
    }
    if(is.null(whichComplex)) return(list(an=an, anno_colors=anno_colors))
    if(is.null(an)) return(NULL)

    anno_colors <- anno_colors[intersect(names(anno_colors),colnames(an))]

    if(length(anno_colors)==0){
        an <- HeatmapAnnotation(df=an, show_legend=show_legend, na_col="white",
                                which=whichComplex,
                                show_annotation_name=show_annotation_name )
    }else{
        an <- HeatmapAnnotation(df=an, show_legend=show_legend, na_col="white",
                                which=whichComplex, col=anno_colors,
                                show_annotation_name=show_annotation_name )
    }
    an
}

.prepData <- function( se, genes=NULL, do.scale=FALSE,
                       assayName=.getDef("assayName"), includeMissing=FALSE ){
    genes <- unique(genes)
    x <- as.matrix(.chooseAssay(se, assayName))
    if(!is.null(genes)) x <- x[intersect(genes,row.names(x)),]
    if(do.scale){
        x <- x[apply(x,1,FUN=sd)>0,]
        x <- t(scale(t(x)))
    }
    if(includeMissing && length(missg <- setdiff(genes, row.names(x)))>0){
        x2 <- matrix( NA_real_, ncol=ncol(x), nrow=length(missg),
                      dimnames=list(missg, colnames(x)) )
        x <- rbind(x,x2)[genes,]
    }
    as.matrix(x)
}

.parseToporder <- function(x, toporder=NULL){
    if(is(x, "SummarizedExperiment")) x <- rowData(x)
    if(is.null(toporder)) return(NULL)
    if(length(toporder)==1 && is.character(toporder)){
        if(toporder %in% colnames(x)){
            toporder <- x[[toporder]]
            names(toporder) <- row.names(x)
        }else{
            stop("Could not interpret `toporder`.")
        }
    }
    if(!is.null(names(toporder))){
        toporder <- toporder[row.names(x)]
    }else{
        names(toporder) <- row.names(x)
    }
    return(toporder)
}

.prepScale <- function(x, hmcols=NULL, breaks=.getDef("breaks")){
    hmcols <- .getHMcols(cols=hmcols)
    if(!is.null(breaks) && !is.na(breaks) && length(breaks)==1 &&
       (!is.logical(breaks) || breaks))
        breaks <- getBreaks(x, length(hmcols)+1, split.prop=breaks)
    if(is.null(breaks) || is.na(breaks) || (is.logical(breaks) && !breaks))
        breaks <- getBreaks(x, length(hmcols)+1, 1, FALSE)
    list(breaks=breaks, hmcols=hmcols)
}

.rbind_all <- function(dfs){
    aac <- unique(unlist(lapply(dfs,colnames)))
    dfs <- lapply(dfs, FUN=function(x){
        x <- as.data.frame(x)
        for(f in setdiff(aac, colnames(x))) x[[f]] <- NA
        x[,aac,drop=FALSE]
    })
    do.call(rbind, dfs)
}
