
#------------------------------------------------
# The following commands are needed to ensure that the roxygen2 package, which deals with documenting the package, does not conflict with the Rcpp package.

# Rcpp          - allows C++ integration
# parallel      - running jobs in parallel
# coda          - MCMC type objects and methods
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
#' @import parallel
#' @import coda
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

#------------------------------------------------
#' @title Bind data to project
#'   
#' @description Load data into a \code{rgeoprofile_project} prior to analysis.
#'   Data must be formatted as a dataframe with samples in rows and loci in
#'   columns. If individuals are polyploid then multiple rows can be used per
#'   sample. Ploidy is allowed to vary between samples, and can be specified in
#'   multiple ways.
#'   
#' @param project an \code{rgeoprofile_project}, as produced by the function
#'   \code{rgeoprofile_project()}
#' @param df a dataframe containing spatial data
#' @param missing_data which value represents missing data
#' @param name optional name of the data set to aid in record keeping
#' @param check_delete_output whether to prompt the user before overwriting 
#'   existing data
#'
#' @export

bind_data <- function(project, df, missing_data = -9, name = NULL, check_delete_output = TRUE) {
  
  # check inputs
  assert_custom_class(project, "rgeoprofile_project")
  assert_dataframe(df)
  
  # check before overwriting existing output
  if (project$active_set>0 && check_delete_output) {
    
    # ask before overwriting. On abort, return original project
    if (!user_yes_no("All existing output and parameter sets for this project will be lost. Continue? (Y/N): ")) {
      return(project)
    }
    
    # replace old project with fresh empty version
    project <- rgeoprofile_project()
  }
  
  # process and perform checks on data
  #dat_processed <- process_data(df, ID_col, pop_col, ploidy_col, data_cols, ID, pop, ploidy, missing_data, wide_format)
  #dat_processed$name <- name
  
  # add data to project
  #project$data <- df
  #project$data_processed <- dat_processed
  
  return(project)
}

