# compute_lethScores.R

#' Computes lethality scores per replicate in an rscreen.data object
#'
#' \code{get.leth.scores} computes lethality scores per replicate in the \code{data_rscreen}
#' object. The \code{data_rscreen} object  may either involve data either for a single screen, or
#' for multiple screens. Note that lethality scores are computed per plate, if a plate variable
#' is available in the 'rscreen.object' given as input. This corrects automatically for any plate
#' effect that affects all wells in a plate to the same extent.

#' @aliases getLethScores
#' @concept scores screen siRNA CRISPR CRISPR-Cas9 crispr-cas9 CRISPR/Cas9 crispr/cas9
#'
#' @param robust logical; \code{TRUE} by default, when the median of the controls is used to
#' compute the scores. If \code{FALSE}, the mean of the controls is used instead.
#' @param rescale logical, indicating whether or not the data is to be rescaled. Defaults to \code{FALSE}.
#' If \code{TRUE}, the screen data will be rescaled per replicate using the function indicated in 'scale.fun'.
#' @param scale.fun the name of the function to be used to rescale the data. Currently the
#' only choices are 'asinh' for the hyperbolic-arc sine transformation, and 'log2' for the
#' log2 transformation. Note that the hyperbolic-arc
#' sine is equivalent to the logarithm for intermediate and high values, and to a linear
#' transformation for low values. This transformation allows for zeros, which is a useful property when
#' transforming count data. Ignored if 'rescale' is \code{FALSE}.
#' @inheritParams get.inv.set
#'
#' @return An object of class rscreen.object, with its data_only slot containing lethality scores.
#' @examples
#' # See vignette
#' @export
#' @seealso \code{\link{read.screen.data}} on reading screen data, \code{\link{combine.screens}}
#' to combine data from multiple screens into a single object, and \code{\link{get.rscreenorm}} to
#' normalize the screen data using the invariant-set quantile normalization.
#'

get.leth.scores <- function(data_rscreen, robust = TRUE, rescale = FALSE,
                            scale.fun = c("asinh", "log2")) {
  if (class(data_rscreen) != "rscreen.object")
    stop("The input data \"data_rscreen\" is not an rscreen.object")
  if (!is.logical(robust))
    stop("robust must be either TRUE or FALSE")
  data_only <- data_rscreen[["data_only"]]

  if(rescale) {
    scfun <- match.arg(tolower(scale.fun), choices = c("asinh", "log2"))
    if(scfun == "asinh") {
      data_only <- asinh(data_only)
  } else if(scfun == "log2") {
      data_only <- log2(data_only)
  }
  }
  data_ann  <- data_rscreen[["data_ann"]]
  use_plate <- data_rscreen[["use_plate"]]
  # Note that data_ann is either a matrix if all experiments/samples are produced with the same protocol,
  # if the input data consists of a single screen (with a number of replicates),
  # or a list with one common annotation matrix "data_ann_common", as well as a "wtype_ann" matrix,
  # with one column per rep and column names corresponding to the 'data_only' column names
  # a "plate_ann" matrix with plate annotation, also with column names as in 'data_only', and a
  # vector 'screen' containing the screen indicators
  # It is possible that wtype_ann is not available, in which case hasWtype = FALSE
  # Then lethality scores are not computed, and the resulting object contains the input values
  let_all <- data_only
  # If hasWtype = FALSE, nothing is done (as no well type information is available)
  if( !hasWtype(data_rscreen) ) {
            warning("No well type information is available, so
                    lethality scores are equal to the input data")
  } else {
  # Single screen
  # Same annotation, as matrix, for all data
  if ( !is.null(dim(data_ann)) ){
    if (use_plate){
      all_plates <- levels(factor(data_ann[, "plate"]))
      for (xj in all_plates) {
        my_ann <-           data_ann[ as.character(data_ann[, "plate" ]) == xj, ]
        dplate <- as.matrix(data_only[as.character(data_ann[, "plate" ]) == xj, ])
        let_all[ as.character(data_ann[ , "plate" ]) == xj, ] <-
          comp.scores(my_data = dplate, wtype_ann = my_ann$wtype, is_robust = robust)
      }
      #       res <- get.list.res(let_all = let_all, my_ann = data_ann$wtype)
    } else if(!use_plate) {
      let_all <- comp.scores(my_data = data_only, wtype_ann = data_ann$wtype, is_robust = robust)
      #         res <- get.list.res(let_all = let_all, my_ann = data_ann$wtype)
    }
    # Multiple screens
    # Different annotation per replicate, in matrix that is part of data_ann, a list
  } else if (class(data_ann) == "list") {
    wtype_mat <- data_ann[["wtype_ann"]]
    # for the entire data at once
    if ( sum(use_plate)==0 ){
      let_all <- comp.scores(my_data = data_only, wtype_ann = wtype_mat, is_robust = robust)
      #        res <- get.list.res(let_all = let_all, my_ann = wtype_mat)
    } else {
      # per plate
      data_ann1 <- list(wtype_ann = data_ann[["wtype_ann"]], plate_ann = data_ann[["plate_ann"]])
      let_all <- comp.scores.diff.ann(my_data = data_only, my_ann = data_ann1, is_robust = robust)
      #        res <- get.list.res(let_all = let_all, my_ann = wtype_mat)
    }
  }
  }
  let_rscreen <- list( data_ann = data_ann, data_only = let_all,
                    use_plate = use_plate )
  class(let_rscreen) <- "rscreen.object"
  let_rscreen
}


#' Computes lethality scores for replicates with the same plate/well annotation
#'
#' \code{comp.scores} is an internal function called by get.leth.scores.
#' It computes lethality scores per replicate, for a data matrix where replicates
#' have the same plate/well annotation.
#'
#' @param my_data data matrix containing only data columns, so all columns are of type numeric.
#' It can be a matrix with a single column, but we must be able to apply nrow, ncol and apply on it.
#' @param wtype_ann object containing the well type annotation. It can be a vector, character or factor.
#' In case different columns of \code{my_data} have different well annotation, this should be a matrix
#' with the same dimensions as \code{my_data}.
#' @param is_robust logical; \code{TRUE} by default, when the median of the controls is used to
#' compute the scores. If \code{FALSE}, the mean of the controls is used instead.
#'
#' @return a matrix with the same dimensions as \code{my_data} containing all lethality scores.
#' @export
#' @seealso \code{\link{read.screen.data}} on reading screen data, \code{\link{combine.screens}}
#' to combine data from multiple screens into a single object, and \code{\link{get.rscreenorm}} to
#' normalize the screen data using the invariant-set quantile normalization.
#'

comp.scores <- function(my_data, wtype_ann, is_robust=TRUE){
  if (is.null( dim(wtype_ann) )) wtype_ann <- matrix( rep( as.character(wtype_ann), ncol(my_data) ),
                                                      nrow = length(wtype_ann), ncol = ncol(my_data) )
  med_neg_mat <- med_pos_mat <- NULL
  if (is_robust) {
    for (xj in 1:ncol(my_data)){
     med_neg_mat <- c(med_neg_mat, median(my_data[ wtype_ann[, xj] == "neg", xj ], na.rm = T))
     med_pos_mat <- c(med_pos_mat, median(my_data[ wtype_ann[, xj] == "pos", xj ], na.rm = T))
    }
  } else {
    for (xj in 1:ncol(my_data)){
      med_neg_mat <- c(med_neg_mat, mean(my_data[ wtype_ann[, xj] == "neg", xj ], na.rm = T))
      med_pos_mat <- c(med_pos_mat, mean(my_data[ wtype_ann[, xj] == "pos", xj ], na.rm = T))
    }
  }
  mat_pos <- matrix(med_pos_mat, nrow = nrow(my_data), ncol = ncol(my_data), byrow = T)
  mat_neg <- matrix(med_neg_mat, nrow = nrow(my_data), ncol = ncol(my_data), byrow = T)
  let_res <- (my_data - mat_neg)/(mat_pos - mat_neg)
  let_res
}

#' Computes lethality scores for replicates with different plate/well annotation
#'
#' \code{comp.scores.diff.ann} is an internal function called by get.leth.scores.
#' It computes lethality scores per replicate and plate in a data matrix, for which plate and well
#' annotation may vary across replicates.
#'
#' @param my_ann list with data frames, each with as many columns as those in \code{my_data} and
#' with the same column names, containing well annotation (wtype_ann) and plate (plate_ann).
#' @inheritParams comp.scores
#'
#' @return a matrix with the same dimensions as \code{my_data} containing all lethality scores.
#' @export
#' @seealso \code{\link{read.screen.data}} on reading screen data, \code{\link{combine.screens}}
#' to combine data from multiple screens into a single object, and \code{\link{get.rscreenorm}} to
#' normalize the screen data using the invariant-set quantile normalization.
#'

comp.scores.diff.ann <- function(my_data, my_ann, is_robust = TRUE){
  wtype_ann <- my_ann[["wtype_ann"]]
  plate_ann <- my_ann[["plate_ann"]]
  all_plates <- vector("list", ncol(my_data))
  for(xj in 1:ncol(my_data)) all_plates[[xj]] <- levels(factor(as.character(plate_ann[, xj])))
  let_res <- my_data
  for(xj in 1:ncol(my_data)) {
    for (xk in all_plates[[xj]]) {
      d_plate <-  matrix( my_data[ plate_ann[, xj] == xk, xj ], nrow = sum( plate_ann[, xj] == xk ), ncol = 1 )
      d_wtype <- wtype_ann[ plate_ann[, xj] == xk, xj ]
      let_res[ plate_ann[, xj] == xk, xj ] <- comp.scores( d_plate, d_wtype, is_robust = is_robust )
    }
  }
  let_res
}

#' Selects lethality scores subsets per feature type
#'
#' \code{get.list.res} is an internal function called by functions involved in constructing
#' invariant sets. It selects subsets of lethality scores per replicate and per feature
#' type, namely library (or sample) feature, positive and negative controls.
#'
#' @param data_only data matrix containing only data columns, so all columns are of type numeric.
#' @param wtype_mat data.frame or matrix with the same dimensions as data_only, where column \emph{k}
#' contains well annotation for column \emph{k} in data_only.
#'
#' @return a list consisting of mat_let, a matrix with lethality scores corresponding to library
#' features; let_neg, a list with as many elements as columns in mat_sample, with negative controls'
#' lethality scores; let_pos, a list with as many elements as columns in mat_sample, with positive
#' controls' lethality scores
#'
#' @export
#' @seealso \code{\link{read.screen.data}} on reading screen data, \code{\link{combine.screens}}
#' to combine data from multiple screens into a single object, and \code{\link{get.rscreenorm}} to
#' normalize the screen data using the invariant-set quantile normalization.
#'

get.list.res <- function(data_only, wtype_mat){
  mat_let <- NULL
  let_pos <- let_neg <- vector("list", ncol(data_only))
  names(let_pos) <- names(let_neg) <- colnames(data_only)
  for (xi in 1:ncol(data_only)) {
    mat_let <- cbind( mat_let, data_only[ wtype_mat[, xi] == "sample", xi] )
    let_neg[[xi]] <- data_only[wtype_mat[, xi] == "neg", xi]
    let_pos[[xi]] <- data_only[wtype_mat[, xi] == "pos", xi]
  }
  list(mat_let = mat_let, let_neg = let_neg, let_pos = let_pos)
}

