#' @title Read Targets File
#'
#' @description \code{importTargets} reads a delimited targets file into a
#'   dataframe.
#'
#' @details The targets file is a text file which contains information about the
#' RNA samples, linking their experimental conditions to their count data files.
#' Rows correspond to samples and columns to their respective count data file
#' locations and the experimental conditions that were applied.  By default, it
#' is expected that the targets file will be tab-delimited but, other delimiters
#' can be used if specified using \code{sep}.
#'
#' The targets file must contain at least three columns with the following,
#' column names:
#'
#' \describe{ \item{\strong{filename}}{name of the file containing the count
#' data for the sample} \item{\strong{condition}}{experimental treatment or
#' condition that was applied} \item{\strong{replicate}}{replicate identifier -
#' can be numeric or character} }
#'
#' \code{importTargets} adds a \code{label} column which contains labels
#' generated using the \code{condition} and \code{replicate} values for the
#' sample to use in plots.
#'
#' @param file A \link{character} string giving the name of the file.
#' @param path A \link{character} string giving the directory containing the
#'   file. This may be omitted if the file is in the current working directory.
#' @param sep A \link{character} used to separate the fields.
#'
#' @family target functions
#' @family import functions
#'
#' @importFrom utils read.table
#' @importFrom methods is
#'
#' @return dataframe: containing the target information and labels
#'
#' @export
#'
#' @examples
#' \dontrun{
#'  # Import tab-delimited targets file from working directory
#'  importTargets("targets.txt")
#'
#'  # Import tab-delimited targets file from another directory
#'  importTargets("targets.txt", path="/path/to/targets")
#'
#'  # Import comma-delimited targets file from working directory
#'  importTargets("targets.txt", sep=",")
#' }

importTargets <- function(file, path=NULL, sep="\t")
{
  stopifnot(is.character(sep))
  if (missing(file)) stop("Could not import targets: need to specify targets file.")
  if (!is.null(path)) file <- file.path(path, file)

  targets <- read.table(file, header=TRUE, stringsAsFactors=FALSE, sep=sep, quote="\'", fill=TRUE, comment.char = "")
  colnames(targets) <- tolower(colnames(targets))

  if(!validateTargets(targets)) stop("Could not import targets: missing required columns.")

  targets$replicate <- gsub(" ", "", targets$replicate, fixed = FALSE)
  targets$condition <- gsub(" ", "", targets$condition, fixed = FALSE)
  targets <- targets[order(targets$condition, targets$replicate), ]

  targets$label <- paste(targets$condition,targets$replicate,sep="_")
  targets$label <- tolower(targets$label)
  if (! "label" %in% colnames(targets)) stop("Could not import targets: label column not generated.")

  return(targets)
}


#' @title Validate Targets File
#'
#' @description Checks whether targets file contains the required columns.
#'
#' @details \code{validateTargets} checks that the targets dataframe contains
#'   the three columns below:
#'
#'   \describe{ \item{\strong{filename}}{name of the file containing the count
#'   data for the sample} \item{\strong{condition}}{experimental treatment or
#'   condition that was applied} \item{\strong{replicate}}{replicate identifier
#'   - can be numeric or character} }
#'
#'   If the target dataframe contains the necessary columns,
#'   \code{validateTargets} returns \code{TRUE}.
#'
#' @param targets  A \link{data.frame} containing the target information.  See
#'   \code{\link{importTargets}} for more information.
#'
#' @family target functions
#' @family import functions
#'
#' @return If valid, returns \code{TRUE}.
#'
#' @export
#'
#' @examples
#'
#' targets <- data.frame( 'filename'=c('file1','file2','file3','file4'),
#'                   'condition' = rep(c('A','B'), times=2),
#'                   'replicate' = c(1,2,1,2))
#'
#' validateTargets(targets)

validateTargets <- function (targets)
{
  if (missing(targets)) stop("Could not import targets: need to specify targets dataframe.")
  stopifnot(is.data.frame(targets))
  if (ncol(targets) < 3) stop("Could not import targets: too few columns in targets dataframe.")

  expected_columns<- c("filename","condition","replicate")
  for (coln in expected_columns)
  {
    if (!coln %in% colnames(targets)) stop(paste0("Could not import targets: ", coln, " column does not exist."))
  }
  return(TRUE)
}
