# get_inv_set.R

#' Wrapper where the core set is created
#'
#' This function yields one core set of lethality scores per replicate, for
#' a given rscreen.object. The core sets are used as basis for quantile-normalizing
#' data for independent screens and their replicates.

#' @aliases getInvSet
#' @concept core set normalization quantile
#'
#' @param data_rscreen an object of class 'rscreen.object', either yielded by \code{\link{read.screen.data}}, or
#' by \code{\link{combine.screens}}. Typically the 'data_only' slot of this object
#' contains lethality scores computed by \code{\link{get.leth.scores}}, which already
#' reflect a value relative to negative and positive controls, per replicate.
#' @param my_gamma scalar or numeric vector, with value(s) between 0 and a finite
#' number. Its default value is \code{NULL}. If given, it will be used to construct the
#' core sets by using the relative distance between negative and positive
#' controls' distributions. If a numeric vector is given, it must contain as many
#' entries as columns in the data_only slot of data_rscreen. A value of 1 means that
#' the core set includes all lethality scores up to the ratio between the MAD
#' of negative controls, and the sum of the MADs of negative and positive controls,
#' assuming the robust variability measure MAD is used.
#' @param var_type string indicating the statistic used as variability measure.
#' Possible values are 'mad' when the median absolute deviation (MAD) is used
#' (the default), 'IQR' when the inter-quartile range is used, and "sd" when the
#' standard deviation is used.
#' @param prop scalar between 0 and 1, corresponding to the desired proportion of
#' lethality scores to be included in the core set, per replicate. The default
#' value is 0.95. If \code{my_gamma} is numeric, any value for this argument is ignored.
#' @param set_shape string indicating the shape of the core set. It accepts
#' "left" (the default), when all scores to the left of the set threshold are
#' included, per replicate, and "centre", when scores included in the core
#' set only exclude the largest and smallest ones, at equal frequencies.
#' @return A logical matrix with the same dimensions as the data_only slot of data_rscreen.
#' @examples
#' # See vignette
#' @export
#' @seealso \code{\link{get.rscreenorm}} to normalize data from multiple screens using core
#' sets and \code{\link{get.leth.scores}} to compute scores that make observations in
#' different screens comparable.
#'

get.inv.set <- function(data_rscreen, my_gamma = NULL, var_type = c("mad", "IQR", "sd"),
                        prop = 0.95, set_shape = c("left","centre")){
  if (class(data_rscreen) != "rscreen.object")
    stop("The input data \"data_rscreen\" is not an rscreen.object")
  if (is.null(my_gamma) & !is.numeric(prop)) stop("You must either define \"mygamma\" or give a numeric value to \"prop\" ")
  var_type <- match.arg(var_type, c("mad", "IQR", "sd"))
  set_shape <- match.arg( tolower(set_shape), c("left","centre","center"))
  data_only <- data_rscreen[["data_only"]]
  if (set_shape == "center") set_shape <- "centre"
    if ( sum(length(prop) %in% c(1, ncol(data_only)))==0 ) stop("Argument \"prop\" must be a scalar or a vector with one entry per data column")
  if ( !is.null(my_gamma)) {
    if ( !is.numeric(my_gamma)) stop("If argument \"my_gamma\" is given, it must be numeric")
    if ( !(length(my_gamma) %in% c(1, ncol(data_only))) ) stop("Argument \"my_gamma\" must be a scalar or a vector with one entry per data column")
  }
  if ( !hasWtype(data_rscreen) ) { data_rscreen[["data_ann"]]$wtype <- rep("sample", nrow(data_only))
    warning("Since no well type information was found, all wells are assumed to be library features, with no assay controls")
  }
  my_th <- get.th.inv.set(data_rscreen = data_rscreen, my_gamma = my_gamma, var_type = var_type, prop = prop, set_shape = set_shape)
  inv_set <-  sel.scores.inv.set(data_rscreen = data_rscreen, th_obj = my_th)
  inv_set
}

#' Extract well type from data_ann and format it as matrix if need be
#'
#' This is an internal function that
#' extracts well type info from 'data_ann' slot. This function is
#' called by \code{\link{get.th.inv.set}}, \code{\link{sel.scores.inv.set}} and others.

#' @aliases getWtype
#' @concept core set normalization quantile threshold
#'
#' @param data_only  data matrix containing lethality scores, per replicate.
#' @param data_ann matrix or list

#' @return A list with the following elements: the first with the thresholds defining
#' the core sets (one value or pair per sample, so either a vector or a matrix),
# and the second the set_shape (a string)
#' @export
#' @seealso \code{\link{get.inv.set}} the wrapper that calls this functions and others
#' to yield the resulting core sets, \code{\link{get.rscreenorm}} to normalize data
#' from multiple screens using core sets and \code{\link{get.leth.scores}} to compute
#' scores that make observations in different screens comparable.
#'

get.wtype <- function(data_only, data_ann)
{
  # Same annotation, as matrix, for all data
  if ( !is.null(dim(data_ann)) ){
    wtype_mat <- matrix(rep(as.character(data_ann$wtype), ncol(data_only)),
                        nrow=nrow(data_only), ncol=ncol(data_only))
    # Different annotation per replicate, in matrix that is part of data_ann, a list
  } else if (class(data_ann) == "list") {
    wtype_mat <- data_ann[["wtype_ann"]]
  }
}

#' Compute thresholds that define core sets per replicate
#'
#' This is an internal function to compute thresholds
#' to be used later on for constructing
#' core sets. Two methods are available: using the distance between negative
#' and positive controls to determine the core set, and including in the set
#' a proportion of the lethality scores. If no well annotation is available indicating
#' which wells correspond to assay controls, then core
#' sets are computed using a proportion of all available features.
#' To be used internally, it is called by
#' get.inv.set.

#' @aliases getThInvSet
#' @concept core set normalization quantile threshold
#'
#' @inheritParams get.inv.set

#' @return A list with the following elements: the first with the thresholds defining
#' the core sets (one value or pair per sample, so either a vector or a matrix),
# and the second the set_shape (a string)
#' @export
#' @seealso \code{\link{get.inv.set}} the wrapper that calls this functions and others
#' to yield the resulting core sets, \code{\link{get.rscreenorm}} to normalize data
#' from multiple screens using core sets and \code{\link{get.leth.scores}} to compute
#' scores that make observations in different screens comparable.
#'

get.th.inv.set <- function(data_rscreen, my_gamma = NULL, var_type = c("mad", "IQR", "sd"),
                        prop = 0.95, set_shape = c("left","centre")){
  data_only <- data_rscreen[["data_only"]]
  set_shape <- match.arg(tolower(set_shape), c("left","centre"))
  if ( !is.null(my_gamma)) {
    if ( !is.numeric(my_gamma)) stop("If argument \"my_gamma\" is given, it must be numeric")
    if ( !(length(my_gamma) %in% c(1, ncol(data_only))) ) stop("Argument \"my_gamma\" must be a scalar or a vector with one entry per data column")
  }
  # Extract wtype_ann: either a matrix or single colum - if the latter, expanded to be a matrix
  # Result: wtype_mat has the same dimensions as data_only
  data_ann  <- data_rscreen[["data_ann"]]
  if ( hasWtype(data_rscreen) ) {
    wtype_mat <- get.wtype(data_only, data_ann)
  } else {
    wtype_mat <- matrix("sample", nrow=nr.screen(data_rscreen), ncol=nc.screen(data_rscreen))
  }
  let_res <- get.list.res(data_only = data_rscreen[["data_only"]], wtype_mat = wtype_mat)
  if (  !is.null(my_gamma) ) {
  var_type <- match.arg(var_type, c("mad", "IQR", "sd"))
  # Compute threshold using gamma, defined as a relation between variability between negative and positive controls
    if (var_type == "mad"){
      mad_neg <- unlist( lapply( let_res[["let_neg"]], mad, na.rm = T ) )
      mad_pos <- unlist( lapply( let_res[["let_pos"]], mad, na.rm = T ) )
    } else if (var_type == "IQR"){
      mad_neg <- unlist( lapply( let_res[["let_neg"]], IQR, na.rm = T ) )
      mad_pos <- unlist( lapply( let_res[["let_pos"]], IQR, na.rm = T ) )
    } else if (var_type == "sd"){
      mad_neg <- unlist( lapply( let_res[["let_neg"]], sd, na.rm = T ) )
      mad_pos <- unlist( lapply( let_res[["let_pos"]], sd, na.rm = T ) )
    }
    # Now to compute the threshold, per screen
      my_th <- mad_neg/( mad_neg + ( mad_pos/my_gamma ) )
  } else if ( set_shape == "left" ){
    # Compute threshold to include proportion 'prop' of lethality scores in core set
    if (length(prop) == 1)  { my_th <- apply(let_res[["mat_let"]], 2, quantile, probs = prop, na.rm = T) }
    if (length(prop) == ncol(data_only)) {
      my_th <- NULL
      for(xj in 1:ncol(data_only)) my_th <- c( my_th, quantile(let_res[["mat_let"]][, xj], probs = prop[xj], na.rm = T) )
    }
  }
    else if ( set_shape == "centre" ) {
      if (length(prop) == 1)  {
        my_th1 <- apply(let_res[["mat_let"]], 2, quantile, probs = (1-prop)/2, na.rm = T)
        my_th2 <- apply(let_res[["mat_let"]], 2, quantile, probs = prop+(1-prop)/2, na.rm = T)
        my_th <- cbind(my_th1, my_th2)
        }
      if (length(prop) == ncol(data_only)) {
        my_th <- NULL
        for(xj in 1:ncol(data_only)) { ####
          my_th1 <- quantile(let_res[["mat_let"]][, xj], probs = (1-prop[xj])/2, na.rm = T)
          my_th2 <- quantile(let_res[["mat_let"]][, xj], probs = prop[xj]+(1-prop[xj])/2, na.rm = T)
          my_th <- cbind(my_th1, my_th2)
          }
        }
    }
  # Threshold and set_shape in list as result
  res <- list( threshold = my_th, set_shape = set_shape )
  res
}




#' Creates an object indicating observations in core set
#'
#' This is an internal function that
#' uses data_rscreen and a threshold object and returns a
#' logical matrix, with the same dimensions as data_only, with \code{TRUE} for every
#' lethality score in the core set, and \code{FALSE} for all others. To be used
#' internally, it is called by get.inv.set.

#' @aliases selScoresInvSet
#' @concept core set normalization quantile threshold
#'
#' @param th_obj object containing thresholds per replicate, as returned by
#' get.th.inv.set
#' @inheritParams get.inv.set
#'
#' @return A logical matrix with the same number of columns as the data_only slot of data_rscreen,
#' and as many rows as annotated with "sample" in the data_only slot
#' @export
#' @seealso \code{\link{get.inv.set}} the wrapper that calls this functions and others
#' to yield the resulting core sets, \code{\link{get.th.inv.set}} that produces the
#' thresholds per replicate, \code{\link{get.rscreenorm}} to normalize data from multiple
#' screens using core sets and \code{\link{get.leth.scores}} to compute scores that
#' make observations in different screens comparable.

sel.scores.inv.set <- function(data_rscreen, th_obj){
  if ( (!is.list(th_obj)) | (length(th_obj) != 2) ) stop("Argument \"th_obj\" should be a list with two elements")
  if ( class(data_rscreen) != "rscreen.object" ) stop("The input data \"data_rscreen\" is not an rscreen.object")
  my_th <- th_obj[["threshold"]]
  set_shape <- th_obj[["set_shape"]]
  data_only <- data_rscreen[["data_only"]]
  nreps <- ncol(data_only)
  # Extract wtype_ann: either a matrix or single colum - if the latter, expanded to be a matrix
  # Result: wtype_ann has the same dimensions as data_only
  data_ann  <- data_rscreen[["data_ann"]]
  if( !hasWtype(data_rscreen) ) {
    wtype_mat <- matrix(rep("sample", nr.screen(data_rscreen)*nreps),
                        nrow = nr.screen(data_rscreen),
                        ncol = nreps)
  } else { wtype_mat <- get.wtype(data_only, data_ann) }
  let_res <- get.list.res(data_only, wtype_mat)
  inv_set <- matrix(FALSE, nrow = nrow(let_res[["mat_let"]]), ncol = nreps)
  colnames(inv_set) <- colnames(data_only)
  if (set_shape == "left") {
    for(xj in 1:nreps)
      inv_set[, xj] <- let_res[["mat_let"]][, xj] <= my_th[xj]
  } else if (set_shape == "centre"){
    for(xj in 1:nreps)
      inv_set[, xj] <- (let_res[["mat_let"]][, xj] >= my_th[xj,1]) & (let_res[["mat_let"]][, xj] <= my_th[xj,2])
    }
 inv_set
}


#' Computes proportion of scores included in core sets
#'
#' This function uses data_rscreen and a threshold object and returns a
#' vector of proportions of scores included in core sets, per replicate. It can also
#' create a barplot of the computed proportions.
#' It is typically used when a \code{my_gamma} is chosen to select core sets,
#' yielding a check on the built sets.

#' @aliases getInvSetProp
#' @concept core set normalization quantile threshold
#'
#' @param th_obj object containing thresholds per replicate, as returned by
#' \code{\link{get.th.inv.set}}.
#' @param plot logical, indicating if a barplot is to be made of the result. Defaults
#' to \code{FALSE}, when no plot is produced.
#' @param cols_list vector of colours to be used for the barplot, so if given it must
#' be as long as the number of columns in data_only slot of data_rscreen. If \code{\link{get.pdata}}
#' was used to read the phenotypic data table, its slot cols_list can be used as input
#' for this argument. By default, a vector of colours created by \code{\link{rainbow}} will
#' be used, a different one for each column in data_only.
#' @inheritParams get.inv.set
#'
#' @return A vector as long as the number of columns in the data_only slot of data_rscreen,
#' with each entry corresponding to the proportion of scores included in an core set.
#' If plot = \code{TRUE}, makes also a barplot of the computed proportions.
#' @examples
#' # See vignette
#' @export
#' @seealso \code{\link{get.inv.set}} the wrapper function that yields core sets,
#' \code{\link{get.th.inv.set}} that produces the thresholds per replicate, \code{\link{get.rscreenorm}}
#' to normalize data from multiple screens using core sets and \code{\link{get.leth.scores}}
#' to compute scores that make observations in different screens comparable.


get.inv.set.prop <- function(data_rscreen, th_obj, plot = FALSE, cols_list = NULL){
  if ( (!is.list(th_obj)) | (length(th_obj) != 2) ) stop("Argument \"th_obj\" should be a list with two elements")
  if ( class(data_rscreen) != "rscreen.object" ) stop("The input data \"data_rscreen\"
                                                is not an rscreen.object")
  nreps <- nc.screen(data_rscreen)
  inv_set <- sel.scores.inv.set(data_rscreen, th_obj)
  prop_inc <- apply(inv_set, 2, mean, na.rm = T)
  if (plot) {
    if (is.null(cols_list))    {
      cols_list <- rainbow( nreps ) } else {
        if ( length(cols_list) < nreps ) {
          cols_list <- rainbow( nreps )
          warning("The length of argument \"cols_list\" is not equal
                   to the total number of replicates, so default colours were used for the plot.")
          }
      }
    barplot(prop_inc, col = cols_list)
  }
  prop_inc
}



