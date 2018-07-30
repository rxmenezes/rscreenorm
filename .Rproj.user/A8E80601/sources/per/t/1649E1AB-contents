# plot_functions.R

#' Creates a colour variable corresponding to a factor.
#'
#' Creates a variable with each entry assuming a different colour,
#' corresponding to levels of a factor. Such a variable is convenient for use with plotting.
#' The colours chosen are picked by \code{\link{rainbow}}, given the number of factor levels.
#'
#' @param myclass a factor, or a variable that can be considered as (and transformed into) a factor.
#' Observations that are NA for the input variable are also NA for the resulting colour variable.
#' @param mystart the starting value for the colour range from \code{\link{rainbow}}. Must be a value between 0 and 1.
#' @param myend the end value for the colour range from \code{\link{rainbow}}. Must be a value between 0 and 1.
#' @return A list with the following slots:  col_var, containing the variable with the colours, with as many entries as
#' the input variable myclass; list_cols, a vector of colours used, as long as the number of levels/distinct values
#' in myclass; and levels_class, a vector containing the categories of myclass, in the same order as in list_cols.
#'
#' @export
#' @seealso \code{\link{rainbow}}
#' @examples
#' # Create a factor
#' fcl <- factor(rep(5:3,each=3))
#' # Now create a vector of colours, one corresponding to each factor level
#' # as well as objects to be used in a legend
#' col.cl <- var.in.colour(fcl)
#' # Make a simple plot displaying the colours and add a legend
#' plot(as.numeric(as.character(fcl)), col=col.cl[["col_var"]], pch=20, cex=2)
#' legend("topright", legend=col.cl[["levels_class"]], pch=20, col=col.cl[["list_cols"]] )


var.in.colour <- function(myclass, mystart = 0.1, myend = 0.9) {
    if (!(is.factor(myclass))) {
        my_fclass <- factor(myclass)
    } else my_fclass <- myclass
    list_cols <- rainbow(nlevels(my_fclass), start = mystart, end = myend)
    my_class_col <- rep(0, length(my_fclass))
    for (xi in 1:length(my_fclass)) {
        if (!is.na(my_fclass[xi])) {
            my_class_col[xi] <- list_cols[my_fclass[xi] == levels(my_fclass)]
        } else {
            my_class_col[xi] <- NA
        }
    }
    list(col_var = my_class_col, list_cols = list_cols, levels_class = levels(my_fclass))
}


#' Plot density of each column in a dataset.
#'
#' Computes a kernel density per column in a dataset and makes a single plot
#' with appropriate limits of all densities together. The graph also includes a legend by default.
#'
#' @param mydata data.matrix which columns are to be made density plots of.
#' @param n  the number of knots at which the density is to be computed - see \code{\link{density}} for details.
#' @param mycols  colours to be assigned to the lines of the densities computed. If left empty, \code{\link{rainbow}} is
#' used to produce one colour per column in mydata.
#' @param na.rm logical - should NAs per column be removed? Its value is passed on to the call to \code{\link{density}}.
#' @param mytitle string containing the title to be used.
#' @param mytype string indicating the line type to be used. If a single value is given, it is recycled across all lines
#' (one for each column in mydata). If a vector is given, it is expected to be as long as the number of columns in mydata.
#' The default is to use solid lines for all densities. For a list of all possible line types, see graphical parameter lty
#' in the help file for \code{\link{par}}.
#' @param myxlab string containing the legend to be used for the x-axis.
#' @param add.leg logical. If \code{TRUE}, adds a legend to the plot.
#' @param myleg  the legend text to be shown. If left empty, and add.leg is\code{TRUE}, the column names of mydata are used.
#' @param xlim a numerical vector of length 2 giving the limits to be used for the x axis.
#' If left empty, the range of values in mydata is used.
#' @return A plot including one density for each column of values in mydata, using plot limits that allow for a complete
#' displ)ay of each density, by default.
#'
#' @export
#' @seealso \code{\link{density}}, \code{\link{rainbow}}
#' @examples
#' mydata <- cbind(replicate(5, rnorm(500)), replicate(5, rnorm(500, mean=1)))
#' colnames(mydata) <- paste("Data", 1:10)
#' pdensity.column(mydata)

pdensity.column <- function(mydata, n=512, mycols=rainbow(ncol(mydata), start=0.1, end=0.9),
                            na.rm=TRUE, mytitle="Data density", mytype = "solid",
                            myxlab="observed values", add.leg = TRUE, myleg=colnames(mydata), xlim=NULL)
{
  if( (length(mytype) > 1) & (length(mytype) != ncol(mydata))
  ) stop("The number of line types must be the same as the number of columns in mydata")
  if(length(mytype) == 1) mytype <- rep(mytype, ncol(mydata))
  if(add.leg & (length(myleg) !=  ncol(mydata))
     ) stop("The number of legend labels given must be equal to the number of columns in mydata")
  if( (!is.null(xlim)) & (length(xlim) != 2)  ) stop("Please provide two values for argument xlim")
  if( (!is.null(xlim)) & (sum(is.na(xlim)) > 0) ) stop("NAs are not allowed in argument xlim")
  myx <- myy <- matrix(0,nrow=n,ncol=ncol(mydata))
  for(xk in 1:ncol(mydata))
  {
    dens.val <- density(mydata[,xk], n=n, na.rm=na.rm)
    myx[,xk] <- dens.val$x
    myy[,xk] <- dens.val$y
  }

  xk <- 1
  if(is.null(xlim)){
    plot(myx[, xk], myy[, xk], type= "l", lty = mytype[xk], col=mycols[xk],
         xlim=range(myx, na.rm= na.rm), ylim=range(myy),
         xlab=myxlab, ylab="density", main=mytitle)
  } else {
    plot(myx[,xk], myy[,xk], type= "l", lty = mytype[xk], col=mycols[xk], xlim=xlim,
         ylim=range(myy), xlab=myxlab, ylab="density", main=mytitle)
  }
  for(xk in 2:ncol(mydata)) lines(myx[,xk], myy[,xk], col=mycols[xk], lty = mytype[xk])
  if(add.leg) legend("topright", legend=myleg, lty=mytype, col=mycols, cex=0.7)
}


#' Plot density of each column in a dataset, using only one well type.
#'
#' Computes a kernel density per column in a dataset and makes a single plot
#' with appropriate limits of all densities together. The graph also includes a legend by default.
#'
#' @param data_rscreen typically an rscreen.object, from which the columns of the 'data_only' slot
#' are to be made density plots of. Otherwise  a numeric matrix or a data frame with only numeric columns,
#' in which case density plots are made of their columns. In the latter case, all entries (rows) are
#' assumed to be of type 'sample'.
#' @param wtype string indicating the well type to be used.
#' Must be a label that exists in the wtype_ann slot.
#' @param n  the number of knots at which the density is to be computed - see \code{\link{density}} for details.
#' @param mycols  colours to be assigned to the lines of the densities computed. If left empty, \code{\link{rainbow}} is
#' used to produce one colour per column in mydata.
#' @param na.rm logical - should NAs per column be removed? Its value is passed on to the call to \code{\link{density}}.
#' @param mytitle string containing the title to be used.
#' @param sel.reps either a character vector containing labels for replicates to be used, or indices of columns from
#' the 'data_only' slot of the data_rscreen object to be used. The default is to use all replicates in the 'data_rscreen' object.
#' @param rescale logical, indicating whether or not the data is to be rescaled. Defaults to \code{FALSE}.
#' If \code{TRUE}, the screen data will be rescaled per replicate using the function indicated in 'scale.fun'.
#' @param scale.fun the name of the function to be used to rescale the data. Currently the
#' only choices are 'asinh' for the hyperbolic-arc sine transformation, and 'log2' for the
#' log2 transformation. Note that the hyperbolic-arc
#' sine is equivalent to the logarithm for intermediate and high values, and to a linear
#' transformation for low values. This transformation allows for zeros, which is a useful property when
#' transforming count data. Ignored if 'rescale' is \code{FALSE}.
#' @param mytype string indicating the line type to be used. If a single value is given, it is recycled across all lines
#' (one for each column in mydata). If a vector is given, it is expected to be as long as the number of columns in mydata.
#' The default is to use solid lines for all densities. For a list of all possible line types, see graphical parameter lty
#' in the help file for \code{\link{par}}.
#' @param myxlab string containing the legend to be used for the x-axis.
#' @param add.leg logical. If \code{TRUE}, adds a legend to the plot.
#' @param myleg  the legend text to be shown. If left empty, and add.leg is \code{TRUE}, the column names of mydata are used.
#' @param leg.side character, the side of the plot where the legend is to be placed. Defaults to 'right',
#' with the other option being 'left'.
#' @param lcex cex parameter for legend, if used. Must be numeric.
#' @param xlim a numerical vector of length 2 giving the limits to be used for the x axis.
#' If left empty, the range of values in mydata is used.
#' @return A plot including one density for each column of values in data_rscreen, using plot limits that allow for a complete
#' displ)ay of each density, by default.
#'
#' @export
#' @seealso \code{\link{density}}, \code{\link{rainbow}}
#' @examples
#' mydata <- cbind(replicate(5, rnorm(500)), replicate(5, rnorm(500, mean=1)))
#' colnames(mydata) <- paste("Data", 1:10)
#' data.ann <- data.frame(geneID = paste("G", 1:500), wtype = rep("sample", 500))
#' data.screen <- list(data_ann = data.ann, data_only = mydata)
#' class(data.screen) <- "rscreen.object"
#' pdw(data.screen, lcex = .8)

pdw <- function(data_rscreen, wtype = "sample", n=512, mycols = NULL,
                na.rm = TRUE, mytitle = paste("Data density", wtype, "wells"),
                sel.reps = scr.names(data_rscreen), rescale = FALSE,
                scale.fun = c("asinh", "log2"), mytype = "solid",
                myxlab = "observed values", add.leg = TRUE,
                myleg = NULL, leg.side = c("right", "left"), lcex = 1, xlim = NULL)
{
  if( class(data_rscreen) != "rscreen.object" ){
    data_only <- as.matrix(data_rscreen)
    data_ann <- data.frame(wtype = rep("sample", nrow(data_only)))
    data_rscreen <- list(data_ann = data_ann,
                      data_only = data_only,
                      use_plate = FALSE)
  }
  leg.side <- match.arg(tolower(leg.side), c("right", "left"))
  if(leg.side == "right") {
    leg1 <- "topright"
  } else {
    leg1 <- "topleft"
  }
  if(is.numeric(sel.reps)) sel.reps <- scr.names(data_rscreen)[sel.reps]
  data_only <- data_rscreen[["data_only"]][, sel.reps, drop = FALSE]
  if(rescale) {
    scfun <- match.arg(tolower(scale.fun), choices = c("asinh", "log2"))
    if(scfun == "asinh") {
      data_only <- asinh(data_only)
    } else if(scfun == "log2") {
      data_only <- log2(data_only)
    }
  }
  data_ann <- data_rscreen[["data_ann"]]
  if(is.null(mycols)) mycols <- rainbow(ncol(data_only), start=0.1, end=0.9)
  if(is.factor(mycols)) mycols <- as.character(mycols)
  if(is.null(myleg)) myleg <- colnames(data_only)
  if( hasWtype(data_rscreen) ){
  if(is.null(dim(data_ann))) { wtype_ann <- data_ann[["wtype_ann"]]
  } else { wtype_ann <- matrix( rep(data_rscreen[["data_ann"]][, "wtype"],
                                ncol(data_only) ),
                                nrow = nrow(data_only),
                                ncol = ncol(data_only) ) }
  } else {
    wtype_ann <- matrix("sample", nrow=nr.screen(data_rscreen), ncol=nc.screen(data_rscreen))
  }
  if( any( apply(wtype_ann == wtype, 2, sum) == 0) ) stop("The wtype given does not exist in at least one replicate")
  xi <- 1
  mydata <-  data_only[ wtype_ann[, xi] == wtype, xi ]
  for(xi in 2:ncol(data_only)) mydata  <- cbind(mydata, data_only[ wtype_ann[, xi] == wtype, xi ])
  if( (length(mytype) > 1) & (length(mytype) != ncol(mydata))
  ) stop("The number of line types must be the same as the number of columns in mydata")
  if(length(mytype) == 1) mytype <- rep(mytype, ncol(mydata))
  if(add.leg & (length(myleg) !=  ncol(mydata))
  ) stop("The number of legend labels given must be equal to the number of columns in mydata")
  if( (!is.null(xlim)) & (length(xlim) != 2)  ) stop("Please provide two values for argument xlim")
  if( !is.null(xlim) ) if( sum(is.na(xlim)) > 0)  stop("NAs are not allowed in argument xlim")
  myx <- myy <- matrix(0,nrow=n,ncol=ncol(mydata))
  for(xk in 1:ncol(mydata))
  {
    dens.val <- density(mydata[,xk], n=n, na.rm=na.rm)
    myx[,xk] <- dens.val$x
    myy[,xk] <- dens.val$y
  }
  xk <- 1
  if(is.null(xlim)){
    plot(myx[, xk], myy[, xk], type= "l", lty = mytype[xk], col=mycols[xk],
         xlim=range(myx, na.rm= na.rm), ylim=range(myy),
         xlab=myxlab, ylab="density", main=mytitle)
  } else {
    plot(myx[,xk], myy[,xk], type= "l", lty = mytype[xk], col=mycols[xk], xlim=xlim,
         ylim=range(myy), xlab=myxlab, ylab="density", main=mytitle)
  }
  for(xk in 2:ncol(mydata)) lines(myx[,xk], myy[,xk], col=mycols[xk], lty = mytype[xk])
  if(add.leg) legend(leg1, legend=myleg, lty=mytype, col=mycols, cex=lcex)
}


#' Selects densities to be used.
#'
#' To be used internally, this function selects columns of a density matrix.
#'
#' @param xi the index of the column to be selected
#' @param mydens a matrix of which columns are to be selected
#' @return The density selected
#'
#' @export
#' @seealso \code{\link{pdmw}}
sel.rep <- function(mydens, xi)
{
  sel.dens <- mydens[, xi]
  sel.dens
}



#' Plot density of each replicate in a dataset, together for multiple well types.
#'
#' Computes a kernel density per replicate and well type in a dataset and makes a single plot
#' with appropriate limits of all densities together, per replicate.
#'
#' @param pdata a pdata object.
#' @param use.wt  character vector giving the unique well type labels to be used to group features.
#' Must contain at least one label, and all labels given
#' must exist in all replicates. The first label is assumed to correspond to library features, internally coded as 'sample'.
#' One density is computed per label in 'use.wt' and all densities for the unique well type labels are
#' plotted in a single graph, for each replicate.
#' @param plot logical, indicating whether or not a plot is to be produced.
#' @param myxlab the label to be used for the x-axis.
#' @param mytitle string containing the title to be used. Its default value is \code{NULL}, in which case the replicate
#' label is used, the same as the column name in the original data. It may also be \code{FALSE}, in which no title is plotted.
#' @param add.text logical, indicating whether a text label is to be added to the plot. If \code{TRUE}, the replicate label will
#' be displayed.
#' @param add.leg logical, indicating whether or not a legend indicating lines corresponding to the well types
#' in \code{use.wt} is to be added.
#' @param col.types character vector, indicating the colours to be used for the second and third well types given in
#' \code{use.wt}. The colour used for the first well type is the one given in the 'cols_list' slot of 'pdata'.
#' @inheritParams pdw
#' @return A plot including one density for each column of values in mydata, using plot limits that allow for a complete
#' displ)ay of each density, by default.
#'
#' @export
#' @seealso \code{\link{density}} for how densities are computed, \code{\link{pdw}} for making plots using values for
#' a single well type per replicate, for all replicates in a single graph.
#' @examples
#' mydata <- cbind(replicate(5, rnorm(500)), replicate(5, rnorm(500, mean=0.5, sd = 2)))
#' mycont <- rbind(replicate(10, rnorm(100)), replicate(10, rnorm(100, mean=1)))
#' mydata <- rbind(mydata, mycont)
#' colnames(mydata) <- paste("Data", 1:10)
#' data.ann <- data.frame(geneID = paste("G", 1:700),
#'                        wtype = c(rep("sample", 500), rep("neg", 100), rep("pos", 100)) )
#' data.screen <- list(data_ann = data.ann, data_only = mydata)
#' class(data.screen) <- "rscreen.object"
#' pdmw(data.screen)

pdmw <- function(data_rscreen, pdata = NULL, use.wt = c("sample", "pos", "neg"), na.rm = TRUE,
                 sel.reps = scr.names(data_rscreen),
                 plot = TRUE, rescale = FALSE, scale.fun = c("asinh", "log2"), myxlab = "",
                 mytitle = NULL, add.text = FALSE, add.leg = TRUE, lcex = 1, col.types = c("red", "blue") ){
  if( class(data_rscreen) != "rscreen.object" ){
    data_only <- as.matrix(data_rscreen)
    data_ann <- data.frame(wtype = rep("sample", nrow(data_only)))
    data_rscreen <- list(data_ann = data_ann,
                      data_only = data_only,
                      use_plate = FALSE)
  }
if(is.numeric(sel.reps)) sel.reps <- colnames(data_rscreen[["data_only"]])[sel.reps]
data_only <- data_rscreen[["data_only"]][, sel.reps, drop = FALSE]
if(rescale) {
  scfun <- match.arg(tolower(scale.fun), choices = c("asinh", "log2"))
  if(scfun == "asinh") {
    data_only <- asinh(data_only)
  } else if(scfun == "log2") {
    data_only <- log2(data_only)
  }
}
if( is.null(pdata) ) {
  pdata <- data.frame( names = colnames(data_only),
                       cols_vec = rainbow(ncol(data_only)) )
  pdata <- list(pdata = pdata,
                cols_list = unique(pdata$cols_vec),
                list_vars = colnames(pdata))
}
mypdata <- list(pdata = pdata[["pdata"]][pdata[["pdata"]]$names %in% sel.reps, ],
                cols_list = pdata[["cols_list"]][pdata[["pdata"]]$names %in% sel.reps],
                list_vars = pdata[["list_vars"]])
data_ann <- data_rscreen[["data_ann"]]
if( hasWtype(data_rscreen) ){
  if(is.null(dim(data_ann))) { wtype_ann <- data_ann[["wtype_ann"]]
  } else { wtype_ann <- matrix( rep(data_rscreen[["data_ann"]][, "wtype"],
                                    ncol(data_only) ),
                                nrow = nrow(data_only),
                                ncol = ncol(data_only) ) }
} else {
  wtype_ann <- matrix("sample", nrow=nr.screen(data_rscreen), ncol=nc.screen(data_rscreen))
}
wsum <- NULL
for(xi in 1:length(use.wt)) wsum <- c( wsum, apply(wtype_ann == use.wt[xi], 2, sum) )
if( any( wsum == 0) ) stop("One or more well type labels given does not exist in at least one replicate")

d.sir.x <- d.sir.y <- NULL
for(xj in 1:ncol(data_only))
{
  dsample <- density(data_only[ wtype_ann[, xj] == use.wt[1], xj], na.rm=na.rm)
  d.sir.x <- cbind(d.sir.x, dsample$x)
  d.sir.y <- cbind(d.sir.y, dsample$y)
}
colnames(d.sir.x) <- colnames(d.sir.y) <- colnames(data_only)
if(length(use.wt)>1){
  d.cont.x <- d.cont.y <- vector("list", length = length(use.wt)-1)
  for(xi in 2:length(use.wt))
  {
    d.conti.x <- d.conti.y <- NULL
    for(xj in 1:ncol(data_only))
    {
      dcont <- density(data_only[ wtype_ann[, xj] == use.wt[xi], xj], na.rm=na.rm)
      d.conti.x <- cbind(d.conti.x, dcont$x)
      d.conti.y <- cbind(d.conti.y, dcont$y)
    }
    colnames(d.conti.x) <- colnames(d.conti.y) <- colnames(data_only)
    d.cont.x[[xi - 1]] <- d.conti.x
    d.cont.y[[xi - 1]] <- d.conti.y
  }
}
all.dens <- vector("list", length = length(use.wt)*2)
all.dens[[1]] <- d.sir.x
all.dens[[2]] <- d.sir.y
for(xi in 2:length(use.wt))
{
  all.dens[[xi*2-1]] <- d.cont.x[[xi-1]]
  all.dens[[xi*2]]   <- d.cont.y[[xi-1]]
}

if(plot){
  for(xi in 1:ncol(data_only))
  {
    data.plot <- lapply(all.dens, sel.rep, xi = xi)
    myxlim <- range(data.plot[[ 1 ]])
    myylim <- range(data.plot[[ 2 ]])
    if(length(use.wt)>1){
      for(xj in 2:length(use.wt))
      {
        myxlim <- range( c( myxlim, range(data.plot[[ xj*2-1 ]]) ) )
        myylim <- range( c( myylim, range(data.plot[[ xj*2 ]]) ) )
      }
    }
    if(is.null(mytitle)) {
      mytitle1 <- mypdata[["pdata"]]$names[xi]
      } else { if(is.logical(mytitle)) { if(!mytitle) { mytitle1 <- ""} } }
    plot(data.plot[[ 1 ]], data.plot[[ 2 ]], type="l", col = cols(mypdata)[xi],
         xlim=myxlim, ylim=myylim, main=mytitle1, # or col=col.cl[levels(f.cl)==xi] ?
         xlab= myxlab, ylab="", lwd=2)
    if(add.text) text(myxlim[1]+0.4*diff(myxlim),
                      myylim[1]+0.95*diff(myylim),
                      labels=mypdata[["pdata"]]$names[xi], pos=4)
    if(length(use.wt)>1){
      for(xj in 2:length(use.wt))
        lines(data.plot[[ xj*2-1 ]], data.plot[[ xj*2 ]], col=col.types[xj-1], lty="dotted")
    }
    if(add.leg) {
      if(length(use.wt)==1) {
        legend("topright", legend=use.wt, lty="solid", lwd=2,
               col = cols(mypdata)[xi], cex=.7)
      } else {
        legend("topright", legend=use.wt, lty=c("solid", rep("dotted",length(use.wt)-1)), lwd=c(2,rep(1, length(use.wt)-1)),
               col=c(cols(mypdata)[xi], col.types[1:(length(use.wt)-1)]), cex=lcex)
      }
    }
  }
}
}
