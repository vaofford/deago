#' @title Write Plot To File
#'
#' @description \code{\link{plotToFile}} writes plot to PNG file.
#'
#' @details If \code{keep_images} is set to 1 in the configuration file (see
#'   \code{\link{importConfig}}) then a subdirectory called \code{images} will
#'   be created in the timestamped results directory (see
#'   \code{\link{makeResultDir}}). \code{\link{plotToFile}} writes the
#'   designated plot to a PNG file within the images subdirectory.
#'
#'   The dimensions of the plot, in pixels, can be defined using \code{width}
#'   and \code{height}.
#'
#' @param plot Plot to be written to file
#' @param resultsDir A \link{character} string giving the path to timestamped
#'   results directory
#' @param file A \link{character} string giving name of file where the plot will
#'   be saved
#' @param width An \link{integer} giving width of the device in pixels
#' @param height An \link{integer} giving height of the device in pixels
#'
#' @importFrom methods is
#' @import ggplot2
#' @importFrom grDevices png dev.off
#'
#' @export

plotToFile <- function(plot, resultsDir, file, width=800, height=600)
{
  if (missing(plot)) stop("Could not write plot: no plot given.")
  if (missing(resultsDir)) stop("Could not write plot: no results directory given.")
  if (missing(file)) stop("Could not write plot: no filename given.")
  stopifnot(is.numeric(width))
  stopifnot(is.numeric(height))
  stopifnot(is.character(resultsDir))
  stopifnot(is.character(file))

  imageDir <- file.path(resultsDir,"images")
  if (dir.exists(imageDir)){
    plotFile = file.path(imageDir, file)
    png(plotFile, width=width, height=height)
    print(plot)
    #ggsave(filename=plotFile, width=width, height=height, units='in')
    dev.off()

    if (!file.exists(plotFile)) warning(paste0("Could not save plot to file: plot file not written: ", plotFile))
  } else {
    stop(paste0("Could not find image directory in ", resultsDir))
  }
}
