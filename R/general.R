#' @title Round Factor Ten
#'
#' @description \code{\link{roundToFactorOfTen}} rounds value to nearest factor of 10.
#'
#' @param x An \link{integer} to be rounded
#'
#' @return An \link{integer} rounded to the nearest factor of ten
#'
#' @export
#'
#' @examples
#' roundToFactorOfTen(19) # would give 20
#' roundToFactorOfTen(194) # would give 200
#' roundToFactorOfTen(1094) # would give 2000

roundToFactorOfTen <- function(x) {
  n=c(1, 2, 4, 5, 6, 8, 10)
  if(length(x) != 1) stop("'x' must be of length 1")
  10^floor(log10(x)) * n[[which(x <= 10^floor(log10(x)) * n)[[1]]]]
}

#' @title Round To Base
#'
#' @description \code{\link{roundToBase}} rounds value to nearest base.
#'
#' @param x An \link{integer} to be rounded
#' @param base An \link{integer} base value
#' @param direction A \link{character} value specifying the rounding directions
#'   (can be either 'up' or 'down')
#'
#' @return An \link{integer} rounded to the nearest base value
#'
#' @export
#'
#' @examples
#' roundToBase(17, 2, 'up') # would give 18
#' roundToBase(17, 10, 'up') # would give 20
#' roundToBase(17, 10, 'down') # would give 10
#' roundToBase(-17, 10, 'down') # would give -20

roundToBase <- function(x, base, direction) {
  if(direction == "up")
  {
    x <- base * ceiling(max(x) / base)
  } else if (direction == "down") {
    x <- base * floor(min(x) / base)
  } else {
    stop("Could not roundToBase, direction not up or down.")
  }
  return(x)
}
