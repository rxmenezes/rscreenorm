# read_screen_data.R
#
#' Reads screen data from a tab-delimited file
#'
#' \code{read.screen.data} reads in a tab-delimited file containing annotation
#' in the first k columns and screen data per replicate in the remaining columns.
#' Annotation columns include gene (target) id as well as feature/well type, where library
#' features, positive and negative controls are identified. Arrayed screens
#' typically also include a plate identifier. All data columns are assumed to have
#' the same annotation. If that is not the case, the data is assumed to correspond to
#' different screens, and it must be loaded for each separate screen first before combining
#' it into a single object using \code{\link{combine.screens}}.

#' @aliases readScreenData
#' @concept input screen siRNA CRISPR CRISPR-Cas9 crispr-cas9 CRISPR/Cas9 crispr/cas9
#'
#' @param filename string, the name of the tab-delimited file containing the data
#' to be read, including the file extension.
#' @param mydir string, indicating the complete folder path where the data file
#' is to be found. It is assumed it ends without any (forward or backward) slashes, for example
#' '/home/projects' instead of '/home/projects/'. If \code{NULL}, it will get the current
#' work directory via \code{\link{getwd}}.
#' @param sep the separator to be used while reading in the data file,
#'        containing a value accepted by read.delim.
#' @param check_names logical, indicating whether or not the column names in the data
#'        file need to be checked for consistency or not. To be used while reading
#'        in the data file. Defaults to \code{FALSE}, when comparability between
#'       column names in the data and in the phenotypic data table is guaranteed.
#' @param ann_cols_names vector of strings containing the names of the annotation
#' columns. It is assumed that these columns are the first few ones in the data file.
#' @param n_ann_cols an integer giving the number of annotation columns (from the
#' first one) in the data read. Any value given is only used if ann_cols_names is
#' \code{NULL}, otherwise it is ignored. If ann_cols_names is given and n_ann_cols
#' is left empty, n_ann_cols is by default the length of ann_cols_names.
#' @param wtype_col string with the name of the annnotation column corresponding
#' to the well annotation column. If not present, the function issues a warning.
#' For lethality scores to be computed, this column must exist.
#' @param wtype_lab labels in wtype_col variable corresponding to library/sample
#' features, positive and negative controls,
#' respectively. Defaults to 'c("sample","pos", "neg")'.
#' @param use_plate logical indicating whether or not the data contains a plate variable.
#' Defaults to \code{FALSE}. If \code{TRUE}, the plate will be used when computing the lethality scores
#' later on with \code{\link{get.leth.scores}}.
#' @param plate_var string with the name of the plate variable.
#'   Only needed if use_plate is \code{TRUE}. Defaults to "plate".
#' @return An object of class rscreen.object.
#' @examples
#' # See vignette
#' @export
#' @seealso \code{\link{combine.screens}} to combine data from multiple screens,
#' \code{\link{combine.pdata}} to combine pdata from multiple screens into a single object,
#' and \code{\link{get.leth.scores}} to compute scores that make observations in
#' different replicates and screens comparable.
#'

read.screen.data <- function(filename, mydir = NULL, sep = "\t",
                             check_names = FALSE, ann_cols_names = NULL,
                             n_ann_cols = length(ann_cols_names),
                             wtype_col = NULL, wtype_lab = c("sample", "pos", "neg"),
                             use_plate = FALSE, plate_var = "plate") {
    if(is.null(mydir)) mydir <- getwd()
    if ((is.null(ann_cols_names)) & (n_ann_cols < 1))
        stop("Please provide either column names for annotation or a non-zero integer corresponding to the number of columns to be used")
    if (!is.logical(check_names)) stop("Argument \"check_names\" must be either TRUE or FALSE")
    data_read <- read.delim(file.path(mydir, filename), sep = sep, check.names = check_names)
    if (length(ann_cols_names) != sum(colnames(data_read) %in% ann_cols_names))
        stop("Column names given in ann_cols_names do not match those in the data file read")
    if ((sum(colnames(data_read) %in% wtype_col)) == 0)
        warning("Column name given by wtype_col does not match any column name in the data file read")
    if (length(wtype_lab) != 3)
        warning("wtype_lab must contain exactly three feature labels")
    if (!is.logical(use_plate)) stop("Argument \"use_plate\" must be either TRUE or FALSE")
    if (!is.null(ann_cols_names)) {
        data_ann <- data_read[, ann_cols_names]
        n_ann_cols <- length(ann_cols_names)
    } else {
        data_ann <- data_read[, 1:n_ann_cols]
    }
#    data_ann[, wtype_col] <- NULL # this messes up the next part, when wtype_col is not NULL
    if(!is.null(wtype_col)) {
    wtype <- as.character(data_ann[, wtype_col])
    wtype[wtype == wtype_lab[1]] <- "sample"
    wtype[wtype == wtype_lab[2]] <- "pos"
    wtype[wtype == wtype_lab[3]] <- "neg"
    data_ann$wtype <- wtype
    hasWtype <- TRUE
    } else { hasWtype <- FALSE }
    if (use_plate){
      if (!is.character(plate_var)) stop("Argument \"plate_var\" must be character")
      if ( sum(colnames(data_ann) %in% plate_var ) == 0 ) stop("Error: \"use_plate\" is TRUE but there is no column called \"use_plate\" amongst the annotation columns")
      data_ann$plate <- data_ann[, plate_var]
      data_ann[, plate_var] <- NULL
    }
    data_only <- as.matrix(data_read[, (n_ann_cols + 1):ncol(data_read)])
    if (!is.numeric(data_only))
        stop("All non-annotation columns in the data read must be numeric")
    data_rscreen <- list(data_ann = data_ann,
                      data_only = data_only,
                      use_plate = use_plate)
    class(data_rscreen) <- "rscreen.object"
    data_rscreen
}



#' Combines different screen data objects into a single one.
#'
#' \code{combine.screens} combines data from multiple screens, already read in
#' using read.screen.data into a single rscreen.object. The different screens
#' involve the same set of library features, but may use different plate sets and/or
#' different (numbers of) controls. As such, different screens share annotation
#' relating to the library features, such as (target) gene id, but are likely to
#' have different plate and/or well annotation. The total number of features
#' measured per screen must be the same for all experiments.
#'
#' @param list_rscreen a list of rscreen objects, such as the ones produced by reading
#' data with read.screen.data.
#' @return An object of class rscreen.object.
#' @examples
#' # See vignette
#' @export
#' @seealso \code{\link{read.screen.data}} to read data for individual screens,
#' \code{\link{combine.pdata}} for combining phenotypic data tables corresponding
#' to different screens that are to be normalized together, and
#' \code{\link{get.leth.scores}} to compute lethality scores for data corresponding
#' to multiple screens
#'
##
# combine.screens
# The output includes a data_ann slot with the common annotation, obtained from the first screen, as well as separate wtype_ann and
# plate_ann matrices with replicate_specific feature annotation
# Note: the function does not check if all annotation columns (except for plate, well)
# are indeed the same across experiments

combine.screens <- function(list_rscreen){
  if (!is.list(list_rscreen)) stop("The argument \"list_rscreen\" must be a list.")
  class_objs <- unlist( lapply(list_rscreen, class) )
  if (any( class_objs != "rscreen.object")) stop("All elements of list_rscreen must be of class rscreen.object")
  n_screens <- length(list_rscreen)
  my_use_plate <- NULL
  for(xi in 1:n_screens) my_use_plate <- c(my_use_plate, list_rscreen[[xi]][["use_plate"]])
  if (length(my_use_plate) < n_screens) stop("Not all screens have valid \"use_plate\" argument")
  if ( (sum(my_use_plate) < n_screens) & (sum(my_use_plate)>0) ) stop("Argument \"use_plate\"must have
                                                                  the same value for all screens")
  xi <- 1
  data_only <- list_rscreen[[ xi ]][["data_only"]]
  nreps <- ncol(data_only)
  screen <- rep(1, nreps)
  data_ann  <- list_rscreen[[ xi ]][["data_ann"]]
  wtype_ann <- matrix(rep(as.character(data_ann[,"wtype"]), ncol(data_only)),
                      nrow=nrow(data_only), ncol=ncol(data_only))
  colnames(wtype_ann) <- colnames(data_only)
  if ( my_use_plate[xi] ){
    plate_ann <- matrix(rep(as.character(data_ann[,"plate"]), ncol(data_only)),
                        nrow=nrow(data_only), ncol=ncol(data_only))
    colnames(plate_ann) <- colnames(data_only)
  }
  data_ann_common <- data_ann[, colnames(data_ann) != "wtype"]
  for(xi in 2:n_screens)
  {
    if (nrow(data_only) != nrow(list_rscreen[[xi]][["data_only"]]))
         stop(paste("Screen",xi," has different total number of features") )
    data_only2 <- list_rscreen[[xi]][["data_only"]]
    nreps2 <- ncol(data_only2)
    screen <- c( screen, rep(xi, nreps2) )
    data_only <- cbind(data_only, data_only2)
    data_ann <- list_rscreen[[xi]][["data_ann"]]
    wtype_ann1 <- matrix(rep(as.character(data_ann[, "wtype"]), ncol(data_only2)),
                        nrow=nrow(data_only2), ncol=ncol(data_only2))
    colnames(wtype_ann1) <- colnames(data_only2)
    wtype_ann <- cbind(wtype_ann, wtype_ann1)
    if ( my_use_plate[xi] ){
      plate_ann1 <- matrix(rep(as.character(data_ann[,"plate"]), ncol(data_only2)),
                          nrow=nrow(data_only2), ncol=ncol(data_only2))
      colnames(plate_ann1) <- colnames(data_only2)
      plate_ann <- cbind(plate_ann, plate_ann1)
    }
    nreps <- nreps+nreps2
  }
  nwt <- apply(apply(wtype_ann, 2, factor), 2, nlevels)
  # The following statement makes sure plate_ann always exists,
  # even when use_plate==FALSE for all screens used
  if(sum(my_use_plate) == 0) plate_ann <- NULL
  data_ann_common$wtype <- data_ann_common$plate <-NULL
  mydata_ann <- list(data_ann = data_ann_common, wtype_ann = wtype_ann, plate_ann = plate_ann, screen = screen)
  res <- list(data_ann = mydata_ann, data_only = data_only, use_plate = my_use_plate)
  class(res) <- "rscreen.object"
  res
}


#' @importFrom grDevices rainbow
#' @importFrom graphics barplot legend lines plot text
#' @importFrom stats IQR density lm mad median quantile sd
#' @importFrom utils read.delim
NULL
