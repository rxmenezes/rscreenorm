# get_pdata.R

#' Reads phenotypic data from a tab-delimited file
#'
#' \code{get.pdata} reads in a tab-delimited file containing sample annotation, such as
#' a phenotypic data table, and creates variables useful for making graphs later on.
#' The input file must have as many rows as replicates included in the study, and as many
#' columns as phenotype variables, which must include at least a variable containing replicate
#' names (the same as given in the data file), one including cell line ids and one including
#' a condition or treatment.

#' @aliases getPdata

#' @param filename string; the name of the tab-delimited file containing the data
#' to be read, including the file extension.
#' @param mydir string; the complete folder path where the data file
#' is to be found. It is assumed it ends without any (forward or backward) slashes, for example
#' '/home/projects' instead of '/home/projects/'. If \code{NULL}, it will get the current
#' work directory via \code{\link{getwd}}.
#' @param sep the separator to be used while reading in the data file,
#'        containing a value accepted by read.delim.
#' @param check_names logical, indicating whether or not the column names in the data
#'        file need to be checked for consistency or not. To be used while reading
#'        in the data file. Defaults to FALSE, when
#'        column names are left unaltered by R.
#' @param names_var string; name of the variable containing the replicate names,
#' corresponding to the same names given in the data file.
#' @param clines_var string; name of the variable containing the cell line names. This may be
#' used later on for graph annotation. This argument is expected to be declared. Of course, the
#' user can instead declare here any variable considered interesting to distinguish sets of
#' replicates, such as the replicate number.
#' @param treat_var string; name of the variable containing the treatment or condition
#' to the well annotation column. Optional.
#' @param sample_ids a character vector containing sample ids yielded by the data file. This
#' is typically the result of colnames of data_only, a slot of the result of read.screen.data. The function
#' checks if the ids given in sample_ids are the same, and are in the same order, as those in names_var.
#' If they are the same but not in the same order, the pdata table has rows re-ordered so that they match
#' the order in sample_ids.  If values in sample_ids are not a permutation of those in names_bar,
#' the function yields an error. We strongly suggest users make use of this variable, but they may ignore
#' it and no check is made.
#' @return A list containing the pdata read with all variables, a colour variable coding for the cell line variable,
#' as well as a vector of the variable names
#' @examples
#' # See vignette
#' @export
#' @seealso \code{\link{read.screen.data}} to read in screen data, \code{\link{combine.pdata}}
#' to combine pdata from multiple screens into a single object, and \code{\link{get.leth.scores}}
#' to compute scores that make observations in different screens comparable.
#'

get.pdata <- function(filename, mydir = NULL, sep = "\t", check_names = FALSE, names_var = "names",
                      clines_var = "clines", treat_var = "treat", sample_ids = NULL) {
  if(is.null(mydir)) mydir <- getwd()
  data_read <- read.delim(file.path(mydir, filename), sep = sep, check.names = check_names)
    if (sum(colnames(data_read) %in% clines_var) == 0)
        stop("The variable name clines_var does not exist in the file read")
    if (sum(colnames(data_read) %in% treat_var) == 0) {
        warning("The variable name treat_var does not exist in the file read")
        treat_var <- NULL
    }
    if (sum(colnames(data_read) %in% names_var) == 0)
        stop("The variable name names_var does not exist in the file read")
    if (!is.null(sample_ids)) {
        if (sum(as.character(data_read[, names_var]) %in% sample_ids) != length(sample_ids))
            stop("The sample_ids given do not all match/are matched by the ids in the names_var variable")
        ids_mat <- data.frame(sample_ids = as.character(sample_ids), orig_order = 1:length(sample_ids), stringsAsFactors = F)
        if (!all(ids_mat$sample_ids == data_read[, names_var])) {
            ids_mat <- ids_mat[order(ids_mat$sample_ids), ]
            data_read <- data_read[order(as.character(data_read[, names_var])), ]
            if (!all(ids_mat$sample_ids == as.character(data_read[, names_var])))
                warning("It was not possible to get sample_ids and names_var in the pheno data in the same order - please check")
            pdata <- data_read[order(ids_mat$orig_order), ]  # Note: pdata is only re-ordered if sample_ids is NULL
        } else {
            pdata <- data_read
        }
    }
    # Make colour vector, with one colour per cell line
    cols_list <- var.in.colour(data_read[, clines_var])
    pdata$cols_vec <- cols_list[[1]]
    list_vars <- c(names_var, clines_var, treat_var)
    # Output object: a list including the pdata, the cols_list list (for legends) and the list of variable names (names/ids, clines, treat_var)
    list(pdata = pdata, cols_list = cols_list, list_vars = list_vars)
}

#' Combines pdata objects corresponding to multiple screens
#'
#' \code{combine.pdata} combines pdata from multiple screens, already read in
#' using get.pdata into a single pdata object. Variables are combined to allow for
#' using them for plots later on. Variable names are standardized to allow for
#' combining all observations smoothly. Other variables that exist in one pdata and not in
#' the other, receive NA for the extra observations.
#'
#' @param list_pdata a list of pdata objects, read in using get.pdata
#' @aliases combinePdata
#' @return A list containing the stacked pdata objects, a single colour variable coding for the cell line variable, as well as
#' a character vector with the variable names containing replicate names, cell line ids and treatment/condition.
#' @examples
#' # See vignette
#' @export
#' @seealso \code{\link{read.screen.data}} to read in screen data, \code{\link{combine.screens}}
#' to combine data from multiple screens,  and \code{\link{get.leth.scores}} to compute scores
#' that make observations in different screens comparable.
#'


combine.pdata <- function(list_pdata){
  n_screens <- length(list_pdata)
  list_vars <- c("names", "clines", "treat")
  xi <- 1
  my_entry <- list_pdata[[xi]]
  my_pdata <- my_entry[["pdata"]]
  my_cols_list <- my_entry[["cols_list"]]
  vars_list <- my_entry[["list_vars"]]
  pdata1 <- my_pdata[, vars_list]
  if ( length(vars_list) < 3)  {
    pdata1$treat <- rep(NA, nrow(pdata1))
  }
  names(pdata1) <- list_vars
  pdata_all <- apply(pdata1, 2, as.character)
if(n_screens >= 2) {
  for(xi in 2:n_screens)
  {
    my_entry <- list_pdata[[xi]]
    my_pdata <- my_entry[["pdata"]]
    my_cols_list <- c(my_cols_list, my_entry[["cols_list"]])
    vars_list <- my_entry[["list_vars"]]
    pdata1 <- my_pdata[, vars_list]
    if ( length(vars_list) < 3)  {
      pdata1$treat <- rep(NA, nrow(pdata1))
    }
    names(pdata1) <- list_vars
    pdata1 <- apply(pdata1, 2, as.character)
    if (ncol(pdata1) < ncol(pdata_all)) stop(paste("pdata", xi, "does not contain all required variables"))
    pdata_all <- rbind(pdata_all, pdata1)
  }
}
  # Require that all colours are different, otherwise make new ones
  if( length(unique(my_cols_list)) < n_screens) {
    cols_list <- var.in.colour(pdata_all[, list_vars[2] ])[[1]]
  } else { cols_list <- my_cols_list }
  pdata_all <- data.frame(pdata_all, check.names=FALSE)
  pdata_all$cols_vec <- cols_list
  # Output object: a list including the pdata, the cols_list list (for legends) and the list of variable names (names/ids, clines, treat_var)
  list(pdata = pdata_all, cols_list = cols_list, list_vars = list_vars)
}


