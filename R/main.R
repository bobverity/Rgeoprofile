
#------------------------------------------------
# The following commands are needed to ensure that the roxygen2 package, which deals with documenting the package, does not conflict with the Rcpp package.

# Rcpp          - allows C++ integration
# fftwtools     - fast Fourier transform, used when smoothing posterior draws into final surface
# ggplot2       - used to produce layered plots
# ggmap         - needed for the get_map function, although ggmap function itself is broken
# RColorBrewer  - used to define default colours in geoPlotAllocation
# rgdal         - required to load shapefiles
# raster        - required when using masks
# viridis       - colour palettes
# ...           - other importFrom declarations recommended by devtools::check

#' @useDynLib RgeoProfile
#' @importFrom Rcpp evalCpp
#' @import fftwtools
#' @import ggplot2
#' @import ggmap
#' @import RColorBrewer
#' @import rgdal
#' @importFrom raster raster extent extent<- rasterize projectRaster distance
#' @import viridis
#' @importFrom grDevices colorRampPalette
#' @import graphics
#' @import stats
#' @import utils
NULL

