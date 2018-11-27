#' @title Make Results Directory
#'
#' @description Creates a timestamped results directory.
#'
#' @details \code{\link{makeResultDir}} creates a timestamped results
#'   directory in the current working directory or in the path specified using
#'   \code{path}.
#'
#'   By default, when the HTML report is generated, the images used in the
#'   report are removed. If \code{keep_images} is set to 1 in the configuration
#'   file (see \code{\link{importConfig}}), then a folder called \code{images}
#'   is created as a subdirectory of the results directory.  All of the images
#'   used in the HTML report will then be stored in this subdirectory.
#'
#' @param path A \code{\link{character}} string giving the location where the
#'   results directory will be created.
#' @param keep_images A \code{\link{logical}} (default=0).
#'
#' @return name of the timestamped results directory
#'
#' @export
#' @examples
#'
#'  # With no image subdirectory
#'  makeResultDir(path=getwd())
#'
#'  # With image subdirectory
#'  makeResultDir(path=getwd(), keep_images=1)
#'

makeResultDir <- function(path=NULL, keep_images=0)
{
  if (is.null(path)) stop("Could not create results directory: must specify output directory.")
  if (!dir.exists(path)) stop(paste0("Could not create results directory: specified path does not exist"))

  resultsDir <- paste0("result_", format( Sys.time(), "%Y%m%d%H%M%S"))
  resultsDir <- file.path(path, resultsDir)
  dir.create(resultsDir)

  checkDirExists(resultsDir)

  if (keep_images == 1)
  {
    imageDir <- file.path(resultsDir,"images")
    dir.create(imageDir)
    checkDirExists(imageDir)
  }

  return(resultsDir)
}

#' @title Check directory exists 
#'
#' @description \code{\link{checkDirExists}} checks expected directory exists.
#'
#' @param dir An \link{character} directory name to be checked
#'
#' @export
#'
#' @examples
#' checkDirExists('dirname')
#' 

checkDirExists <- function(path)
{
  if (!dir.exists(path)) stop(paste0(path," does not exist."))
}
