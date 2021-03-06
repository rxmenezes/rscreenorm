% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_functions.R
\name{pdensity.column}
\alias{pdensity.column}
\title{Plot density of each column in a dataset.}
\usage{
pdensity.column(mydata, n = 512, mycols = rainbow(ncol(mydata), start = 0.1,
  end = 0.9), na.rm = TRUE, mytitle = "Data density", mytype = "solid",
  myxlab = "observed values", add.leg = TRUE, myleg = colnames(mydata),
  xlim = NULL)
}
\arguments{
\item{mydata}{data.matrix which columns are to be made density plots of.}

\item{n}{the number of knots at which the density is to be computed - see \code{\link{density}} for details.}

\item{mycols}{colours to be assigned to the lines of the densities computed. If left empty, \code{\link{rainbow}} is
used to produce one colour per column in mydata.}

\item{na.rm}{logical - should NAs per column be removed? Its value is passed on to the call to \code{\link{density}}.}

\item{mytitle}{string containing the title to be used.}

\item{mytype}{string indicating the line type to be used. If a single value is given, it is recycled across all lines
(one for each column in mydata). If a vector is given, it is expected to be as long as the number of columns in mydata.
The default is to use solid lines for all densities. For a list of all possible line types, see graphical parameter lty
in the help file for \code{\link{par}}.}

\item{myxlab}{string containing the legend to be used for the x-axis.}

\item{add.leg}{logical. If \code{TRUE}, adds a legend to the plot.}

\item{myleg}{the legend text to be shown. If left empty, and add.leg is\code{TRUE}, the column names of mydata are used.}

\item{xlim}{a numerical vector of length 2 giving the limits to be used for the x axis.
If left empty, the range of values in mydata is used.}
}
\value{
A plot including one density for each column of values in mydata, using plot limits that allow for a complete
displ)ay of each density, by default.
}
\description{
Computes a kernel density per column in a dataset and makes a single plot
with appropriate limits of all densities together. The graph also includes a legend by default.
}
\examples{
mydata <- cbind(replicate(5, rnorm(500)), replicate(5, rnorm(500, mean=1)))
colnames(mydata) <- paste("Data", 1:10)
pdensity.column(mydata)
}
\seealso{
\code{\link{density}}, \code{\link{rainbow}}
}
