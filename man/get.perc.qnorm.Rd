% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_qnorm.R
\name{get.perc.qnorm}
\alias{get.perc.qnorm}
\alias{getpercqorm}
\title{Quantile-normalizes percentiles for observations in the core set}
\usage{
get.perc.qnorm(data_perc)
}
\arguments{
\item{data_perc}{a numeric matrix of percentiles, with as many columns as samples and as many rows
as percentiles}
}
\value{
a numeric matrix with the same dimensions as the input \code{data_perc}, containing the quantile-
normalized percentiles for the core sets, per replicate
}
\description{
This is an internal function that yields
quantile-normalizes percentiles for the data in the core set.
This function is called by both \code{\link{get.rscreenorm}} and \code{\link{get.qnorm}}. It is
to be used internally, so has no functionality other than what is described.
}
\seealso{
\code{\link{get.rscreenorm}} which yields core sets and quantile-normalizes data at
once, code{\link{get.qnorm}} which quantile-normalizes data for a given core set,
\code{\link{get.leth.scores}} to compute lethality scores, and \code{\link{get.perc}} to compute
the percentiles.
}
\concept{
scores screen siRNA CRISPR CRISPR-Cas9 crispr-cas9 CRISPR/Cas9 crispr/cas9 quantile normalization
}
