# Auxiliary functions to be used by many other functions in the package

#' Checks if an rscreen.object has a well type variable
#'
#' To be used internally, this function checks for a well type variable
#' within the rscreen.object
#'
#' @param data_rscreen an rscreen.object
#' @return A logical value that is TRUE if a well type variable is found.
#' @export

hasWtype <- function(data_rscreen){
  if (!is.null(dim(data_rscreen[["data_ann"]]))) { haswtype <- "wtype" %in% colnames(data_rscreen[["data_ann"]])
  } else {
    haswtype <- "wtype_ann" %in% names(data_rscreen[["data_ann"]])
  }
haswtype
}

#' Gets the number of features in the rscreen.object
#'
#' To be used internally, this function extracts the number of rows
#' in the data_only slot of  the rscreen.object
#'
#' @param data_rscreen an rscreen.object
#' @return A numeric value.
#' @export

nr.screen <- function(data_rscreen){
  nrow(data_rscreen[["data_only"]])
}

#' Gets the number of replicate in the rscreen.object
#'
#' To be used internally, this function extracts the number of columns
#' in the data_only slot of  the rscreen.object.
#'
#' @param data_rscreen an rscreen.object
#' @return A numeric value.
#' @export

nc.screen <- function(data_rscreen){
  ncol(data_rscreen[["data_only"]])
}


#' Extracts the number of well types
#'
#' To be used internally, this function extracts the number of levels
#' in the well type variable within a given rscreen.object. If the object
#' has no well type, it will return 0.
#'
#' @param data_rscreen an rscreen.object
#' @return A numeric value.
#' @export

nwt <- function(data_rscreen){
  if(hasWtype(data_rscreen)) {
    if(wtype %in% names(data_rscreen[["data_ann"]])) {
      nwt <- nlevels(factor(data_rscreen[["data_ann"]][, "wtype"]))
    } else {
      nwt <- apply(apply(wtype_ann, 2, factor), 2, nlevels)
    }
  } else { nwt <- 0}
  nwt
}

#' Extracts sample names of an rscreen object
#'
#' This function extracts the column names
#' of the 'data_only' slot of  a given rscreen.object. If 'the argument 'data_rscreen'
#' is not of class 'rscreen.object', the column names of the object are returned.
#'
#' @param data_rscreen an rscreen.object
#' @return A character vector.
#' @export

scr.names <- function(data_rscreen){
  if( class(data_rscreen) != "rscreen.object" ){
    vnames <- colnames(data_rscreen)
  } else{
    vnames <- colnames(data_rscreen[["data_only"]])
  }
  vnames
  }

#' Extracts the screen data of an rscreen object
#'
#' This function extracts the  'data_only' slot of  a given rscreen.object.
#'
#' @param data_rscreen an rscreen.object
#' @return A numeric data matrix.
#' @export

scr.data <- function(data_rscreen){
   data_rscreen[["data_only"]]
}


#' Extracts the annotation data of an rscreen object
#'
#' This function extracts the  'data_only' slot of  a given rscreen.object.
#'
#' @param data_rscreen an rscreen.object
#' @return A data frame containing annotation columns.
#' @export

scr.ann <- function(data_rscreen){
  data_rscreen[["data_ann"]]
}

#' Extracts the vector of colours from a pdata object
#'
#' This function extracts the column 'cols_vec' from the 'pdata slot of
#' a given pdata object.
#'
#' @param pdata a pdata object consisting of a list including a data frame called 'pdata',
#' containing phenotypic information (on columns) for samples (on rows). One of those is a
#' column called 'cols_vec' containing colours that are to be used in plots, corresponding
#' to the samples.
#' @return A character vector containing colours corresponding to samples in the pdata object.
#' @export

cols <- function(pdata){
  as.character(pdata[["pdata"]]$cols_vec)
}



#' Splits and extracts single value from each entry in a vector of strings
#'
#' To be used internally, this function splits a string vector using a given separator,
#' extracts one of the results and combines all extracted results in a vector. The split
#' is done by using \code{\link{strsplit}}. This function
#' was written to be used sequentially on a vector via \code{\link{sapply}}.
#'

#' @param xi  the element to be selected, typically is a number varying from 1 to the length of the character vector.
#' @param char.vector the character vector containing the strings to be split.
#' @param split.by the string to be used to split each entry in char.vector - see \code{\link{strsplit}} for details.
#' @param result.sel integer indicating the index of the element after splitting to be selected.
#' @param fixed logical, to be passed on to \code{\link{strsplit}}.
#' \code{\link{strsplit}} returns  a list of results, of which the one with index result.sel is selected.
#' @return A  single string.
#' @examples
#' xs <- paste(LETTERS[1:5], 5:1, sep="_")
#' xs
#' sapply(1:length(xs), mysplit, char.vector = xs, split.by = "_", result.sel = 1)
#' # Compare it with the result from directly applying strsplit:
#' strsplit(xs, "_")
#' @export

mysplit <- function(xi, char.vector, split.by, result.sel, fixed=FALSE){
  char.vector <- as.character(char.vector)
  mychar <- char.vector[xi]
  # If we want to have all bits of the split as a matrix/array, we can use
  # split.result <- strsplit(mychar,split=split.by)[[result.sel]]
  if(!fixed) {split.result <- strsplit(mychar,split=split.by)[[1]][result.sel]
  split.result} else {
    split.result <- strsplit(mychar,split=split.by,fixed=fixed)[[1]][result.sel]
    split.result
  }
}



