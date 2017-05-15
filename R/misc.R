
# -----------------------------------
# The following commands are needed to ensure that the roxygen2 package, which deals with documenting the package, does not conflict with the Rcpp package.

# TODO - check that packages rgeos and RgoogleMaps are in fact needed. These used to be included in the ringHS function, so might be needed there.

#' @useDynLib RgeoProfile
#' @importFrom Rcpp evalCpp
#' @import fftwtools
#' @import ggplot2
#' @import ggmap
#' @import RColorBrewer
#' @import rgl
#' @import sp
#' @import rgeos
#' @import RgoogleMaps
NULL
