% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/misc.R
\name{getBreaks}
\alias{getBreaks}
\title{getBreaks}
\usage{
getBreaks(x, n, split.prop = 0.98, symmetric = TRUE)
}
\arguments{
\item{x}{A matrix of log2FC (or any numerical values centered around 0)}

\item{n}{The desired number of breaks.}

\item{split.prop}{The proportion of the data points to plot on a linear
scale; the remaining will be plotted on a scale with regular frequency per
step (quantile).}

\item{symmetric}{Logical; whether breaks should be symmetric around 0
(default TRUE)}
}
\value{
A vector of breaks of length = `n`
}
\description{
Produces symmetrical breaks for a color scale, with the scale steps
increasing for large values, which is useful to avoid outliers influencing
too much the color scale.
}
\examples{
dat <- rnorm(100,sd = 10)
getBreaks(dat, 10)
}
