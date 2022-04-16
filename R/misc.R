
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
        if(!is.null(fromAssay) && grepl("^log",fromAssay, ignore.case=TRUE)){
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
	    if(toAssay=="log2FC")
	        assays(x)$scaledLFC <- sechm::safescale(assays(x)$log2FC,
	                                                center=FALSE, byRow=TRUE)
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
#' @param norm Logical; whether to calculate logcpm (TMM normalization).
#' @param lfc_group the colData column to use to calculate foldchange. If
#' NULL (default), no foldchange assay will be computed.
#'
#' @return A SummarizedExperiment
#' @importFrom edgeR cpm calcNormFactors DGEList
#' @import SummarizedExperiment
#' @importFrom S4Vectors metadata metadata<- SimpleList
#' @export
flattenPB <- function(pb, norm=TRUE, lfc_group=NULL){
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
#' @param filename xlsx file name
#' @param addSheets An optional list of additional tables to save as sheets.
#'
#' @return Saves to file.
#'
#' @examples
#' data("SE", package="SEtools")
#' # not run
#' # se2xls(SE, filename="SE.xlsx")
#'
#' @importFrom openxlsx write.xlsx
#' @export
se2xls <- function(se, filename, addSheets=NULL){
    a <- list( sample_annotation=as.data.frame(colData(se)) )
    if(ncol(rowData(se))>0) a$feature_annotation=as.data.frame(rowData(se))
    a <- c(a, as.list(assays(se)), addSheets)
    write.xlsx(a, file=filename, row.names=TRUE, col.names=TRUE)
}


.prepareAnnoDF <- function(an, anno_colors, fields, whichComplex=NULL,
                           show_legend=TRUE, show_annotation_name=TRUE,
                           dropEmptyLevels=TRUE, anno_name_side=NULL){
    if(!is.null(whichComplex))
        whichComplex <- match.arg(whichComplex, c("row","column"))
    an <- an[,intersect(fields, colnames(an)),drop=FALSE]
    an <- as.data.frame(an)
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
            }else if(!is.null(anno_colors[[i]]) &&
                     !is.function(anno_colors[[i]])){
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
                                which=whichComplex, annotation_name_side=anno_name_side,
                                show_annotation_name=show_annotation_name )
    }else{
        an <- HeatmapAnnotation(df=an, show_legend=show_legend, na_col="white",
                                which=whichComplex, col=anno_colors,
                                annotation_name_side=anno_name_side,
                                show_annotation_name=show_annotation_name )
    }
    an
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

