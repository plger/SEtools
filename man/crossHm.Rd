% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/crossHm.R
\name{crossHm}
\alias{crossHm}
\title{crossHm}
\usage{
crossHm(
  ses,
  genes,
  do.scale = TRUE,
  uniqueScale = FALSE,
  assayName = .getDef("assayName"),
  sortBy = seq_along(ses),
  only.common = TRUE,
  cluster_cols = FALSE,
  cluster_rows = is.null(sortBy),
  toporder = NULL,
  hmcols = NULL,
  breaks = .getDef("breaks"),
  gaps_at = .getDef("gaps_at"),
  gaps_row = NULL,
  anno_rows = .getDef("anno_rows"),
  anno_columns = .getDef("anno_columns"),
  name = NULL,
  anno_colors = list(),
  show_rownames = NULL,
  merge_legends = FALSE,
  show_colnames = FALSE,
  rel.width = NULL,
  ...
)
}
\arguments{
\item{ses}{A (named) list of
\code{\link[SummarizedExperiment]{SummarizedExperiment-class}}.}

\item{genes}{A vector of genes/row.names to plot.}

\item{do.scale}{Logical; whether to scale rows in each SE (default TRUE).}

\item{uniqueScale}{Logical; whether to force the same colorscale for
each heatmap.}

\item{assayName}{The name of the assay to use; if multiple names are given,
the first available will be used. Defaults to "logcpm", "lognorm".}

\item{sortBy}{Names or indexes of `ses` to use for sorting rows (default all)}

\item{only.common}{Logical; whether to plot only rows common to all SEs
(default TRUE).}

\item{cluster_cols}{Logical; whether to cluster columns (default FALSE).}

\item{cluster_rows}{Logical; whether to cluster rows (default TRUE if
`do.sortRows=FALSE`, FALSE otherwise).}

\item{toporder}{Optional verctor of categories on which to supra-order when
sorting rows, or name of a `rowData` column to use for this purpose.}

\item{hmcols}{Colors for the heatmap.}

\item{breaks}{Breaks for the heatmap colors. Alternatively, symmetrical
breaks can be generated automatically by setting `breaks` to a numerical
value between 0 and 1. The value is passed as the `split.prop` argument to
the \code{\link{getBreaks}} function, and indicates the proportion of the
points to map to a linear scale, while the more extreme values will be
plotted on a quantile scale. `breaks=FALSE` will disable symmetrical scale
and quantile capping, while retaining automatic breaks. `breaks=1` will
produce a symmetrical scale without quantile capping.}

\item{gaps_at}{Columns of `colData` to use to establish gaps between columns.}

\item{gaps_row}{A named vector according to which rows will be split.}

\item{anno_rows}{Columns of `rowData` to use for annotation.}

\item{anno_columns}{Columns of `colData` to use for annotation.}

\item{name}{The title of the heatmap key.}

\item{anno_colors}{List of colors to use for annotation.}

\item{show_rownames}{Whether to show row names (default TRUE if 50 rows or
less).}

\item{merge_legends}{Logical; passed to
\code{\link[ComplexHeatmap]{draw-HeatmapList-method}}}

\item{show_colnames}{Whether to show column names (default FALSE).}

\item{rel.width}{Relative width of the heatmaps}

\item{...}{Any other parameter passed to each call of
\code{\link[ComplexHeatmap]{Heatmap}}.}
}
\value{
A Heatmap list.
}
\description{
These functions have been moved and will be deprecated from this package;
please use the
[sechm](https://bioconductor.org/packages/release/bioc/html/sechm.html)
package instead.
}
\details{
Plot a multi-panel heatmap from a list of
\code{\link[SummarizedExperiment]{SummarizedExperiment-class}}.
}
\examples{
data("SE", package="SEtools")
se1 <- SE[,1:10]
se2 <- SE[,11:20]
se3 <- mergeSEs( list(se1=se1, se2=se2) )

}
