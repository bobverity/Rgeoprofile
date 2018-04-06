
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
#' @importFrom raster raster
#' @import viridis
#' @importFrom grDevices colorRampPalette
#' @importFrom graphics abline axis box legend par points polygon segments
#' @importFrom stats dnorm optim rnorm
#' @importFrom utils data flush.console
NULL

#------------------------------------------------
#' Create Rgeoprofile data object
#'
#' Simple function that ensures that input data is in the correct format required by Rgeoprofile. Takes longitude and latitude as input vectors and returns these same values in list format. If no values are input then default values are used.
#'
#' @param longitude the locations of the observed data in degrees longitude.
#' @param latitude the locations of the observed data in degrees latitude.
#'
#' @export
#' @examples
#' # John Snow cholera data
#' data(Cholera)
#' geoData(Cholera$longitude, Cholera$latitude)
#' 
#' # simulated data
#' sim <-rDPM(50, priorMean_longitude = -0.04217491, priorMean_latitude = 
#' 51.5235505, alpha=1, sigma=1, tau=3)
#' geoData(sim$longitude, sim $latitude)

geoData <- function(longitude=NULL, latitude=NULL) {
	
	# use example data if none read in
	if (is.null(longitude) & is.null(latitude)) {
    #data(LondonExample_crimes)
    longitude <- LondonExample_crimes$longitude
    latitude <- LondonExample_crimes$latitude
	} else {
		if (is.null(longitude) | is.null(latitude)) {
			stop("Both longitude and latitude arguments must be used, or alternatively both arguments must be set to NULL to use default values")
		}
	}
	
	# combine and return
	ret <- list(longitude=longitude, latitude=latitude)
	return(ret)
}

#------------------------------------------------
#' Create sources data object in same format as observations
#'
#' Simple function that ensures that sources are in the correct format required by Rgeoprofile. Takes longitude and latitude as input vectors and returns these same values in list format. If no values are input then default values are used.
#'
#' @param longitude the locations of the potential sources in degrees longitude.
#' @param latitude the locations of the potential sources in degrees latitude.
#'
#' @export
#' @examples
#' # John Snow cholera data
#' data(WaterPumps)
#' geoDataSource(WaterPumps$longitude, WaterPumps$latitude)
#' 
#' # simulated data
#' sim <-rDPM(50, priorMean_longitude = -0.04217491, priorMean_latitude = 
#' 51.5235505, alpha=1, sigma=1, tau=3)
#' geoDataSource(sim$source_longitude, sim$source_latitude)

geoDataSource <- function(longitude=NULL, latitude=NULL) {
    
  # generate dummy data if none read in
  if (is.null(longitude) & is.null(latitude)) {
		#data(LondonExample_sources)
  	longitude <- LondonExample_sources$longitude
		latitude <- LondonExample_sources$latitude
  } else {
  	if (is.null(longitude) | is.null(latitude))
  		stop("Both longitude and latitude arguments must be used, or alternatively both arguments must be set to NULL to use default values")
  }
  
  # combine and return
  ret <- list(longitude=longitude, latitude=latitude)
  return(ret)
}

#------------------------------------------------
#' Create Rgeoprofile parameters object
#'
#' This function can be used to generate parameters in the format required by other Rgeoprofile functions. Parameter values can be specified as input arguments to this function, or alternatively if data is input as an argument then some parameters can take default values directly from the data.
#'
#' @param data observations in the format defined by geoData().
#' @param sigma_mean the mean of the prior on sigma (sigma = standard deviation of the dispersal distribution).
#' @param sigma_var the variance of the prior on sigma.
#' @param sigma_squared_shape as an alternative to defining the prior mean and variance of sigma, it is possible to directly define the parameters of the inverse-gamma prior on sigma^2. If so, this is the shape parameter of the inverse-gamma prior.
#' @param sigma_squared_rate the rate parameter of the inverse-gamma prior on sigma^2.
#' @param priorMean_longitude the mean longitude of the normal prior on source locations (in degrees). If NULL then defaults to the midpoint of the range of the data, or -0.1277 if no data provided.
#' @param priorMean_latitude the mean latitude of the normal prior on source locations (in degrees). If NULL then defaults to the midpoint of the range of the data, or 51.5074 if no data provided.
#' @param tau the standard deviation of the normal prior on source locations, i.e. how far we expect sources to lie from the centre. If NULL then defaults to the maximum distance of any observation from the prior mean, or 10.0 if no data provided.
#' @param alpha_shape shape parameter of the gamma prior on the parameter alpha.
#' @param alpha_rate rate parameter of the gamma prior on the parameter alpha.
#' @param chains number of MCMC chains to use in the burn-in step.
#' @param burnin number of burn-in iterations to be discarded at start of MCMC.
#' @param samples number of sampling iterations. These iterations are used to generate final posterior distribution.
#' @param burnin_printConsole how frequently (in iterations) to report progress to the console during the burn-in phase.
#' @param samples_printConsole how frequently (in iterations) to report progress to the console during the sampling phase.
#' @param longitude_minMax vector containing minimum and maximum longitude over which to generate geoprofile. If NULL then defaults to the range of the data plus a guard rail on either side, or c(-0.1377,-0.1177) if no data provided.
#' @param latitude_minMax vector containing minimum and maximum latitude over which to generate geoprofile. If NULL then defaults to the range of the data plus a guard rail on either side, or c(51.4974, 51.5174) if no data provided.
#' @param longitude_cells number of cells in the final geoprofile (longitude direction). Higher values generate smoother distributions, but take longer to run.
#' @param latitude_cells number of cells in the final geoprofile (latitude direction). Higher values generate smoother distributions, but take longer to run.
#' @param guardRail when data input is used, longitude_minMax and latitude_minMax default to the range of the data plus a guard rail. This parameter defines the size of the guard rail as a proportion of the range. For example, a value of 0.05 would give an extra 5 percent on the range of the data.
#'
#' @export
#' @examples
#' # John Snow cholera data
#' data(Cholera)
#' d <- geoData(Cholera$longitude, Cholera$latitude)
#' # define parameters such that the model fits sigma from the data
#' geoParams(data = d, sigma_mean = 1.0, sigma_squared_shape = 2, 
#' chains = 10, burnin = 1000, samples = 10000, guardRail = 0.1)
#' 
#' # simulated data
#' sim <-rDPM(50, priorMean_longitude = -0.04217491, priorMean_latitude = 
#' 51.5235505, alpha=1, sigma=1, tau=3)
#' d <- geoData(sim$longitude, sim $latitude)
#' # use a fixed value of sigma
#' geoParams(data = d, sigma_mean = 1.0, sigma_var = 0,
#' chains=10, burnin=1000, samples = 10000, guardRail = 0.1)

geoParams <- function(data=NULL, sigma_mean=1, sigma_var=NULL, sigma_squared_shape=NULL, sigma_squared_rate=NULL, priorMean_longitude=NULL, priorMean_latitude=NULL, tau=NULL, alpha_shape=0.1, alpha_rate=0.1, chains=10, burnin=1e3, samples=1e4, burnin_printConsole=100, samples_printConsole=1000, longitude_minMax=NULL, latitude_minMax=NULL, longitude_cells=500, latitude_cells=500, guardRail=0.05) {
  
  # if data argument used then get prior mean and map limits from data
  if (!is.null(data)) {
    
    # check correct format of data
    geoDataCheck(data, silent=TRUE)
    
    # if prior mean not defined then set as midpoint of data
    if (is.null(priorMean_longitude)) { priorMean_longitude <- sum(range(data$longitude))/2 }
    if (is.null(priorMean_latitude)) { priorMean_latitude <- sum(range(data$latitude))/2 }
    
    # convert data to bearing and great circle distance, and extract maximum great circle distance to any point
    data_trans <- latlon_to_bearing(priorMean_latitude, priorMean_longitude, data$latitude, data$longitude)
    dist_max <- max(data_trans$gc_dist)
    
    # use maximum distance as default value of tau
    if (is.null(tau)) { tau <- dist_max }
    
    if (is.null(longitude_minMax)) {
      lon_range <- diff(range(data$longitude))
      lon_min <- min(data$longitude) - guardRail*lon_range
      lon_max <- max(data$longitude) + guardRail*lon_range
      longitude_minMax <- c(lon_min, lon_max)
    }
    if (is.null(latitude_minMax)) {
      lat_range <- diff(range(data$latitude))
      lat_min <- min(data$latitude) - guardRail*lat_range
      lat_max <- max(data$latitude) + guardRail*lat_range
      latitude_minMax <- c(lat_min, lat_max)
    }
  }
  
  # if data argument not used then set prior values from input arguments, or use defaults if NULL
  if (is.null(data)) {
    if (is.null(priorMean_longitude)) { priorMean_longitude <- -0.1277 }
    if (is.null(priorMean_latitude)) { priorMean_latitude <- 51.5074 }
    if (is.null(longitude_minMax)) { longitude_minMax <- priorMean_longitude + c(-0.01,0.01) }
    if (is.null(latitude_minMax)) { latitude_minMax <- priorMean_latitude + c(-0.01,0.01) }
    if (is.null(tau)) { tau <- 10 }
  }
  
  # initialise shape and rate parameters for prior on sigma^2
  alpha <- NULL
  beta <- NULL
  
  # calculate shape and rate parameters of prior on sigma-squared (ie. alpha and beta) from input arguments. The user has several options:
  #   1) specify sigma_mean and sigma_var, in which case alpha and beta are calculated from these values
  #   2) specify sigma_mean but not sigma_var, in which case alpha must also be specified
  #   3) specify neither sigma_mean nor sigma_var, in which case alpha and beta must be specified
  
  #-------- option 1
  if (!is.null(sigma_mean) & !is.null(sigma_var)) {
    
    # if using fixed sigma model then no need to calculate alpha and beta. Otherwise use values of sigma_mean and sigma_var to search for the unique alpha and beta that define the distribution
    if (sigma_var==0) {
      cat('Using fixed sigma model')
    } else {
      cat('Using sigma_mean and sigma_var to define prior on sigma')
      ab <- get_alpha_beta(sigma_mean, sigma_var)
      alpha <- ab$alpha
      beta <- ab$beta
    }
    
  #-------- option 2
  } else if (!is.null(sigma_mean) & is.null(sigma_var)) {
    
    if (!is.null(sigma_squared_shape)) {
      cat('Using sigma_mean and sigma_squared_shape to define prior on sigma')
      alpha <- sigma_squared_shape
      if (alpha<=1) { stop('sigma_squared_shape must be >1') }
      beta <- exp(2*log(sigma_mean) + 2*lgamma(alpha) - 2*lgamma(alpha-0.5))
      epsilon <- sqrt(beta)*gamma(alpha-0.5)/gamma(alpha)
      sigma_var <- beta/(alpha-1)-epsilon^2
    }
    
  #-------- option 3
  } else if (is.null(sigma_mean) & is.null(sigma_var)) {
    
    if (!is.null(sigma_squared_shape) & !is.null(sigma_squared_rate)) {
      cat('Using sigma_squared_shape and sigma_squared_rate to define prior on sigma')
      alpha <- sigma_squared_shape
      beta <- sigma_squared_rate
      sigma_mean <- sqrt(beta)*gamma(alpha-0.5)/gamma(alpha)
      sigma_var <- beta/(alpha-1)-sigma_mean^2
    }
    
  }
  
  # check that chosen inputs do in fact uniquely define the distribution. At this stage alpha and beta are only allowed to be NULL under the fixed-sigma model
  if (is.null(alpha) | is.null(beta)) {
    returnError <- TRUE
    if (!is.null(sigma_var)) {
      if (sigma_var==0) {
        returnError <- FALSE
      }
    }
    if (returnError) {
      stop("Current prior parameters on sigma do not fully specify the distribution. Must specify either 1) a prior mean and variance on sigma, 2) a prior mean on sigma and a prior shape on sigma^2, 3) a prior shape and prior rate on sigma^2.")
    }
  }
  
  # set model parameters
  model <- list(sigma_mean=sigma_mean, sigma_var=sigma_var, sigma_squared_shape=alpha, sigma_squared_rate=beta, priorMean_longitude=priorMean_longitude, priorMean_latitude=priorMean_latitude, tau=tau, alpha_shape=alpha_shape, alpha_rate=alpha_rate)
  
  # set MCMC parameters
  MCMC <- list(chains=chains, burnin=burnin, samples=samples, burnin_printConsole=burnin_printConsole, samples_printConsole=samples_printConsole)
  
  # set output parameters    
  output <- list(longitude_minMax=longitude_minMax, latitude_minMax=latitude_minMax, longitude_cells=longitude_cells, latitude_cells=latitude_cells)
  
  # combine and return
  ret <- list(model=model, MCMC=MCMC, output=output)
  return(ret)
}

#------------------------------------------------
#' Import shapefile
#'
#' TODO
#'
#' @param fileName TODO
#'
#' @export
#' @examples
#' # TODO

geoShapefile <- function(fileName=NULL) {
  
  # load north London boroughs by default
  if (is.null(fileName)) {
    fileName <- system.file('extdata', 'London_north', package='RgeoProfile')
  }
  
  # load shapefile
  ret <- rgdal::readOGR(fileName, verbose=FALSE)
  
  return(ret)
}

#------------------------------------------------
#' Check data
#'
#' Check that all data for use in Rgeoprofile MCMC is in the correct format.
#'
#' @param data a data list object, as defined by geoData().
#' @param silent whether to report if data passes checks to console.
#'
#' @export
#' @examples
#' # John Snow cholera data
#' data(Cholera)
#' d <- geoData(Cholera$longitude, Cholera$latitude)
#' geoDataCheck(d)
#' 
#' # simulated data
#' sim <-rDPM(50, priorMean_longitude = -0.04217491, priorMean_latitude = 
#' 51.5235505, alpha=1, sigma=1, tau=3)
#' d <- geoData(sim$longitude, sim $latitude)
#' geoDataCheck(d)

geoDataCheck <- function(data, silent=FALSE) {
    
  # check that data is a list
  stopifnot(is.list(data))
  
  # check that contains longitude and latitude
  stopifnot("longitude" %in% names(data))
  stopifnot("latitude" %in% names(data))
  
  # check that data values are correct format and range
  stopifnot(is.numeric(data$longitude))
  stopifnot(all(is.finite(data$longitude)))
  stopifnot(is.numeric(data$latitude))
  stopifnot(all(is.finite(data$latitude)))
  
  # check same number of observations in logitude and latitude, and n>1
  stopifnot(length(data$longitude)==length(data$latitude))
  stopifnot(length(data$longitude)>1)
  
  # if passed all checks
  if (!silent) { cat("data file passed all checks\n") }
}

#------------------------------------------------
#' Check parameters
#'
#' Check that all parameters for use in Rgeoprofile MCMC are in the correct format.
#'
#' @param params a list of parameters, as defined by geoParams().
#' @param silent whether to report passing check to console.
#'
#' @export
#' @examples
#' # John Snow cholera data
#' data(Cholera)
#' d <- geoData(Cholera$longitude, Cholera$latitude)
#' p <- geoParams(data = d, sigma_mean = 1.0, sigma_squared_shape = 2)
#' geoParamsCheck(p)
#' 
#' # simulated data
#' sim <-rDPM(50, priorMean_longitude = -0.04217491, priorMean_latitude = 
#' 51.5235505, alpha=1, sigma=1, tau=3)
#' d <- geoData(Cholera$longitude, Cholera$latitude)
#' p <- geoParams(data = d, sigma_mean = 1.0, sigma_var=0)
#' geoParamsCheck(p)

geoParamsCheck <- function(params, silent=FALSE) {
  
  # check that params is a list
  if (!is.list(params))
    stop("params must be in list format")
      
  # check that contains 'model', 'MCMC' and 'output' as sublists
  if (!"model"%in%names(params) | !"MCMC"%in%names(params) | !"output"%in%names(params))
    stop("params must contain sublists 'model', 'MCMC' and 'output'")

  # check that 'model', 'MCMC' and 'output' are indeed lists
  if (!is.list(params$model))
    stop("params$model must be in list format")
  if (!is.list(params$MCMC))
    stop("params$MCMC must be in list format")
  if (!is.list(params$output))
    stop("params$output must be in list format")

  #---------------------------------------

  # check that params$model contains all necessary parameters
  if (!("sigma_mean"%in%names(params$model)))
    stop("params$model must contain parameter 'sigma_mean'")
  if (!("sigma_var"%in%names(params$model)))
    stop("params$model must contain parameter 'sigma_var'")
  if (!("sigma_squared_shape"%in%names(params$model)))
    stop("params$model must contain parameter 'sigma_squared_shape'")
  if (!("sigma_squared_rate"%in%names(params$model)))
    stop("params$model must contain parameter 'sigma_squared_rate'")
  if (!("priorMean_longitude"%in%names(params$model)))
    stop("params$model must contain parameter 'priorMean_longitude'")
  if (!("priorMean_latitude"%in%names(params$model)))
    stop("params$model must contain parameter 'priorMean_latitude'")
  if (!("tau"%in%names(params$model)))
    stop("params$model must contain parameter 'tau'")
  if (!("alpha_shape"%in%names(params$model)))
    stop("params$model must contain parameter 'alpha_shape'")
  if (!("alpha_rate"%in%names(params$model)))
    stop("params$model must contain parameter 'alpha_rate'")

  # check that params$model values are correct format and range
  if (!is.numeric(params$model$sigma_mean) | !is.finite(params$model$sigma_mean))
    stop("params$model$sigma_mean must be numeric and finite")
  if (params$model$sigma_mean<=0)
    stop("params$model$sigma_mean must be greater than 0")
  if (!is.numeric(params$model$sigma_var) | !is.finite(params$model$sigma_var))
    stop("params$model$sigma_var must be numeric and finite")
  if (params$model$sigma_var<0)
    stop("params$model$sigma_var must be greater than or equal to 0")
  
  # the only time that sigma_squared_shape and sigma_squared_rate are allowed to be NULL is under the fixed sigma model
  if (is.null(params$model$sigma_squared_shape) | is.null(params$model$sigma_squared_rate)) {
	if (params$model$sigma_var!=0) {
		stop('params$model$sigma_squared_shape and params$model$sigma_squared_rate can only be NULL under the fixed sigma model, i.e. when params$model$sigma_var==0. ')
	}
  }

  if (!is.null(params$model$sigma_squared_shape)) {
  	if (!is.numeric(params$model$sigma_squared_shape) | !is.finite(params$model$sigma_squared_shape))
    	stop("params$model$sigma_squared_shape must be numeric and finite")
  }
  if (!is.null(params$model$sigma_squared_rate)) {
  	if (!is.numeric(params$model$sigma_squared_rate) | !is.finite(params$model$sigma_squared_rate))
    	stop("params$model$sigma_squared_rate must be numeric and finite")
  }
  if (!is.numeric(params$model$priorMean_longitude) | !is.finite(params$model$priorMean_longitude))
    stop("params$model$priorMean_longitude must be numeric and finite")
  if (!is.numeric(params$model$priorMean_latitude) | !is.finite(params$model$priorMean_latitude))
    stop("params$model$priorMean_latitude must be numeric and finite")
  if (!is.numeric(params$model$tau) | !is.finite(params$model$tau))
    stop("params$model$tau must be numeric and finite")
  if (params$model$tau<=0)
    stop("params$model$tau must be greater than 0")
  if (!is.numeric(params$model$alpha_shape) | !is.finite(params$model$alpha_shape))
    stop("params$model$alpha_shape must be numeric and finite")
  if (params$model$alpha_shape<=0)
    stop("params$model$alpha_shape must be greater than 0")
  if (!is.numeric(params$model$alpha_rate) | !is.finite(params$model$alpha_rate))
    stop("params$model$alpha_rate must be numeric and finite")
  if (params$model$alpha_rate<=0)
    stop("params$model$alpha_rate must be greater than 0")

  #---------------------------------------

  # check that params$MCMC contains all necessary parameters
  if (!("chains"%in%names(params$MCMC)))
    stop("params$MCMC must contain parameter 'chains'")
  if (!("burnin"%in%names(params$MCMC)))
    stop("params$MCMC must contain parameter 'burnin'")
  if (!("samples"%in%names(params$MCMC)))
    stop("params$MCMC must contain parameter 'samples'")
  if (!("burnin_printConsole"%in%names(params$MCMC)))
    stop("params$MCMC must contain parameter 'burnin_printConsole'")
  if (!("samples_printConsole"%in%names(params$MCMC)))
    stop("params$MCMC must contain parameter 'samples_printConsole'")
  
  # check that params$MCMC values are correct format and range
  if (!is.numeric(params$MCMC$chains) | !is.finite(params$MCMC$chains))
    stop("params$MCMC$chains must be numeric and finite")
  if (params$MCMC$chains<=1)
    stop("params$MCMC$chains must be 2 or more")
  if (!is.numeric(params$MCMC$burnin) | !is.finite(params$MCMC$burnin))
    stop("params$MCMC$burnin must be numeric and finite")
  if (params$MCMC$burnin<0)
    stop("params$MCMC$burnin must be greater than or equal to 0")
  if (!is.numeric(params$MCMC$samples) | !is.finite(params$MCMC$samples))
    stop("params$MCMC$samples must be numeric and finite")
  if (params$MCMC$samples<=0)
    stop("params$MCMC$samples must be greater than 0")
  if (!is.numeric(params$MCMC$burnin_printConsole) | !is.finite(params$MCMC$burnin_printConsole))
    stop("params$MCMC$burnin_printConsole must be numeric and finite")
  if (params$MCMC$burnin_printConsole<=0)
    stop("params$MCMC$burnin_printConsole must be greater than 0")
  if (!is.numeric(params$MCMC$samples_printConsole) | !is.finite(params$MCMC$samples_printConsole))
    stop("params$MCMC$samples_printConsole must be numeric and finite")
  if (params$MCMC$samples_printConsole<=0)
    stop("params$MCMC$samples_printConsole must be greater than 0")
  
  #---------------------------------------
  
  # check that params$output contains all necessary parameters
  if (!("longitude_minMax"%in%names(params$output)))
    stop("params$output must contain parameter 'longitude_minMax'")
  if (!("latitude_minMax"%in%names(params$output)))
    stop("params$output must contain parameter 'latitude_minMax'")
  if (!("longitude_cells"%in%names(params$output)))
    stop("params$output must contain parameter 'longitude_cells'")
  if (!("latitude_cells"%in%names(params$output)))
    stop("params$output must contain parameter 'latitude_cells'")
  
  # TODO: SOME MORE CHECKS ON FORMAT OF params$output?
  
  #---------------------------------------
  
  # if passed all checks
  if (!silent) { cat("params file passed all checks\n") }
}

#------------------------------------------------
#' MCMC under Rgeoprofile model
#'
#' This function carries out the main MCMC under the Rgeoprofile model. Posterior draws are smoothed to produce a posterior surface, and converted into a geoProfile. Outputs include posterior draws of alpha and sigma under the variable-sigma model.
#'
#' @param data input data in the format defined by geoData().
#' @param params input parameters in the format defined by geoParams().
#' @param lambda bandwidth to use in posterior smoothing. If NULL then optimal bandwidth is chosen automatically by maximum-likelihood.
#'
#' @export
#' @examples
#' # John Snow cholera data
#' data(Cholera)
#' d <- geoData(Cholera$longitude, Cholera$latitude)
#' p <- geoParams(data = d, sigma_mean = 1.0, sigma_squared_shape = 2)
#' m <- geoMCMC(data = d, params = p, lambda=0.05)
#' 
#' # simulated data
#' sim <-rDPM(50, priorMean_longitude = -0.04217491, priorMean_latitude = 
#' 51.5235505, alpha=1, sigma=1, tau=3)
#' d <- geoData(sim$longitude, sim $latitude)
#' p <- geoParams(data = d, sigma_mean = 1.0, sigma_squared_shape = 2)
#' m <- geoMCMC(data = d, params = p)

geoMCMC <- function(data, params, lambda=NULL) {
    
    # check that data and parameters in correct format
    geoDataCheck(data)
    geoParamsCheck(params)
    cat("\n")
    
    # extract ranges etc. from params object
    min_lon <- params$output$longitude_minMax[1]
    max_lon <- params$output$longitude_minMax[2]
    min_lat <- params$output$latitude_minMax[1]
    max_lat <- params$output$latitude_minMax[2]
    cells_lon <- params$output$longitude_cells
    cells_lat <- params$output$latitude_cells
    cellSize_lon <- (max_lon-min_lon)/cells_lon
    cellSize_lat <- (max_lat-min_lat)/cells_lat
    breaks_lon <- seq(min_lon, max_lon, l=cells_lon+1)
    breaks_lat <- seq(min_lat, max_lat, l=cells_lat+1)
    mids_lon <- breaks_lon[-1] - cellSize_lon/2
    mids_lat <- breaks_lat[-1] - cellSize_lat/2
    mids_lon_mat <- outer(rep(1,length(mids_lat)), mids_lon)
    mids_lat_mat <- outer(mids_lat, rep(1,length(mids_lon)))
    
    # transform data to cartesian coordinates relative to centre of prior. After transformation data are defined relative to point 0,0 (i.e. the origin represents the centre of the prior). Add transformed coordinates to data object before feeding into C++ function
    data_cartesian <-latlon_to_cartesian(params$model$priorMean_latitude, params$model$priorMean_longitude, data$latitude, data$longitude)
    data$x <- data_cartesian$x
    data$y <- data_cartesian$y
    
    # if using fixed sigma model then change alpha and beta from NULL to -1. This value will be ignored, but needs to be numeric before feeding into the C++ function.
    if (params$model$sigma_var==0) {
        params$model$sigma_squared_shape <- -1
        params$model$sigma_squared_rate <- -1
    }
    
    # carry out MCMC using efficient C++ function
    rawOutput <- C_geoMCMC(data, params)
    
    # extract mu draws and convert from cartesian to lat/lon coordinates
    mu_draws <- cartesian_to_latlon(params$model$priorMean_latitude, params$model$priorMean_longitude, rawOutput$mu_x, rawOutput$mu_y)
    
    # produce smoothed surface
    mu_smooth <- geoSmooth(mu_draws$longitude, mu_draws$latitude, breaks_lon, breaks_lat, lambda)
    
    # calculate coordinates of lat/lon matrix in original cartesian coordinates
    cart <-latlon_to_cartesian(params$model$priorMean_latitude, params$model$priorMean_longitude, mids_lat_mat, mids_lon_mat)
    
    # produce prior matrix. Note that each cell of this matrix contains the probability density at that point multiplied by the size of that cell, meaning the total sum of the matrix from -infinity to +infinity would equal 1. However, as the matrix is limited to the region specified by the limits, in reality this matrix will usually sum to less than 1.
    priorMat <- dnorm(cart$x, sd=params$model$tau) * dnorm(cart$y, sd=params$model$tau) * (cellSize_lon*cellSize_lat)
    
    # combine prior surface with stored posterior surface (the prior never fully goes away under the DPM model)
    n <- length(data$longitude)
    alpha <- rawOutput$alpha
    posteriorMat <-  mu_smooth + priorMat*mean(alpha/(alpha+n))
    
    # produce geoprofile
    gp <- geoProfile(posteriorMat)
    
    # calculate posterior allocation
    allocation <- matrix(unlist(rawOutput$allocation), n, byrow=T)
    allocation <- data.frame(allocation/params$MCMC$samples)
    names(allocation) <- paste("group", 1:ncol(allocation), sep="")
    
    # get single best posterior grouping
    bestGrouping <- apply(allocation, 1, which.max)
    
    # calculate posterior co-allocation
    coAllocation <- matrix(unlist(rawOutput$coAllocation), n, byrow=T)/params$MCMC$samples
    diag(coAllocation) <- 1
    coAllocation[row(coAllocation)>col(coAllocation)] <- NA
    
    # finalise output format
    output <- list()
    output$priorSurface <-  priorMat
    output$posteriorSurface <-  posteriorMat
    output$geoProfile <-  gp
    output$midpoints_longitude <- mids_lon
    output$midpoints_latitude <- mids_lat
    output$sigma <- rawOutput$sigma
    output$alpha <- alpha
    output$allocation <- allocation
    output$bestGrouping <- bestGrouping
    output$coAllocation <- coAllocation
    
    return(output)
}

#------------------------------------------------
#' Calculate geoprofile from surface
#'
#' Converts surface to hitscore percentage
#'
#' @param surface matrix to convert to geoprofile
#'
#' @export
#' @examples
#' # John Snow cholera data
#' data(Cholera)
#' d <- geoData(Cholera$longitude, Cholera$latitude)
#' p <- geoParams(data = d, sigma_mean = 1.0, sigma_squared_shape = 2)
#' m <- geoMCMC(data = d, params = p, lambda=0.05)
#' gp <- geoProfile(m$posteriorSurface)
#' 
#' # simulated data
#' sim <-rDPM(50, priorMean_longitude = -0.04217491, priorMean_latitude = 
#' 51.5235505, alpha=1, sigma=1, tau=3)
#' d <- geoData(sim$longitude, sim $latitude)
#' p <- geoParams(data = d, sigma_mean = 1.0, sigma_squared_shape = 2)
#' m <- geoMCMC(data = d, params = p)
#' gp <- geoProfile(m$posteriorSurface)

geoProfile <- function(surface) {
  
  # check that surface is in correct format
  stopifnot(is.matrix(surface))
  
  # create geoprofile from surface
  ret <- matrix(rank(surface, ties.method="first"), nrow=nrow(surface), byrow=FALSE)
  ret[is.na(surface)] <- NA
  ret <- 100 * (1 - (ret-1)/max(ret, na.rm=TRUE))
  
  return(ret)
}

#------------------------------------------------
#' Calculate hitscores
#'
#' Calculate hitscores of the potential sources for a given surface (usually the geoprofile).
#'
#' @param params TODO
#' @param source longitude and latitude of one or more source locations in the format defined by geoDataSource().
#' @param surface TODO
#'
#' @export
#' @examples
#' # TODO
#' 
geoReportHitscores <- function(params, source, surface) {
  
  # get size of cells
  delta_lat <- diff(params$output$latitude_minMax)/params$output$latitude_cells
  delta_lon <- diff(params$output$longitude_minMax)/params$output$longitude_cells
  
  # get index of closest point to each source
  index_lat <- round((source$latitude - params$output$latitude_minMax[1])/delta_lat + 0.5)
  index_lon <- round((source$longitude - params$output$longitude_minMax[1])/delta_lon + 0.5)
  
  # return hitscore percentages in data frame
  hs <- surface[cbind(index_lat, index_lon)]
  ret <- data.frame(longitude=source$longitude, latitude=source$latitude, hitscorePercentage=hs)
  
  return(ret)
}

#------------------------------------------------
#' Extract latitude and longitude of points identified as sources by geoMCMC()
#' 
#' This function takes the output of geoMCMC() and, for each 'crime', extracts the group to which it is assigned with the highest probability. For each group, the model returns the mean lat/long of all crimes assigned to that group.
#' 
#' @param mcmc Model output in the format produced by geoMCMC().
#' @param data Crime site data, in the format produced by geoData().
#' 
#' @export
#' @examples
#' # simulated data
#' sim <-rDPM(50, priorMean_longitude = -0.04217491, priorMean_latitude = 
#' 51.5235505, alpha=10, sigma=1, tau=3)
#' d <- geoData(sim$longitude, sim $latitude)
#' s <- geoDataSource(sim$source_lon, sim$source_lat)
#' p <- geoParams(data = d, sigma_mean = 1.0, sigma_squared_shape = 2)
#' m <- geoMCMC(data = d, params = p)
#' # extract sources identified by the model
#' ms <- geoModelSources(mcmc = m, data = d)
#' # plot data showing the sources identified by the model (note: NOT the actual suspect sites)
#' geoPlotMap(data = d, source = ms$model_sources, params = p, breakPercent = seq(0, 10, 1), 
#'                   mapType = "roadmap", surfaceCols =c("red", "orange","yellow","white"),
#'                   crimeCol = "black", crimeCex = 2, sourceCol = "red", sourceCex = 2,
#'                   surface = m$geoProfile, gpLegend=TRUE, opacity = 0.4)

geoModelSources <- function (mcmc, data) {
  
  # get mean over data, split by best group
  lon <- mapply(mean, split(data$longitude, mcmc$bestGrouping))
  lat <- mapply(mean, split(data$latitude, mcmc$bestGrouping))
  
  return(list(longitude=lon, latitude=lat))
}

#------------------------------------------------
#' Calculate and plot hit scores based on a ring search
#'
#' Second attempt at ring search without using other packages.
#'
#' @param params Parameters list in the format defined by geoParams().
#' @param data Data object in the format defined by geoData().
#' @param source Potential sources object in the format defined by geoDataSource().
#' @param mcmc mcmc object of the form produced by geoMCMC(). 
#'
#' @export
#' @examples
#' # TODO

geoRing <- function(params, data, source, mcmc) {
  
  # Calculates the percentage of the grid that must be searched before reaching each source under a ring search strategy. This search strategy assumes that we start from a given crime and search outwards in a circle of increasing radius until we reach a source. As there are multiple crimes the strategy assumes a separate individual searching from each crime simultaneously at an equal rate.
  # The basic logic of the approach here is that calculating the final search radius needed to identify a source is easy - it is simply the minimum radius from any crime to this source. The difficulty is calculating the amount of grid that will have been explored by the time we reach this radius, as circles will often overlap and the intersection should not be double-counted (we assume searching moves on if the area has already been explored by someone else). This is done by brute force - a grid is created and cells are filled in solid if they have been explored. The total percentage of filled cells gives the hitscore percentage. The distance matrices used in this brute force step are needed repeatedly, and so they are computed once at the begninning to save time.
  
  # get number of crimes and sources
  n <- length(data$latitude)
  ns <- length(source$source_latitude)
  
  # create matrices giving lat/lon at all points in search grid
  lonVec <- mcmc$midpoints_longitude
  latVec <- mcmc$midpoints_latitude
  lonMat <- matrix(rep(lonVec,each=length(latVec)), length(latVec))
  latMat <- matrix(rep(latVec,length(lonVec)), length(latVec))
  
  # calculate great circle distance from every data point to every point in search grid. This list of distance matrices will be used multiple times so best to pre-compute here.
  ret <- matrix(-Inf, nrow=nrow(lonMat), ncol=ncol(lonMat))
  for (i in 1:n) {
    neg_dist_i <- -latlon_to_bearing(data$latitude[i], data$longitude[i], latMat, lonMat)$gc_dist
    ret[neg_dist_i>ret] <- neg_dist_i[neg_dist_i>ret]
  }
  
  # return negative distance
  return(ret)
}

#------------------------------------------------
#' Incorporate shapefile or raster information into a geoprofile
#' 
#' This function allows information from a shapefile or raster to be incorporated within the geoprofile. For example, we might wish to exclude areas not on land, or weight the probabilities within a specific postcode differently. The spatial object used should be a SpatialPolygonsDataFrame as produced by the package sp or a raster. 
#' 
#' @param probSurface the original geoprofile, usually the object $posteriorSurface produced by geoMCMC().
#' @param params an object produced by geoParams().
#' @param mask the spatial information to include. Must be one of SpatialPolygonsDataFrame, SpatialLinesDataFrame or RasterLayer.
#' @param scaleValue different functions depending on value of "operation". For "inside' or "outside", the value by which probabilities should be multiplied inside or outside the shapefile; set to zero to exclude completely. For "near" and "far", the importance of proximity to, or distance from, the object described in the RasterLayer or SpatialPointsDataFrame. Not used for "continuous".
#' @param operation thow to combine the surface and the new spatial information. Must be one of "inside", "outside", "near", "far" or "continuous". The first two multiply areas inside or outside the area described in the shapefile (or raster) by scaleValue. "near" or "far" weight the geoprofile by its closeness to (or distance from) the area described in the shapefile (or raster). Finally, "continuous" uses a set of numerical values (eg altitude) to weight the geoprofile.
#' @param maths one of "add", "subtract", multiply" or "divide. The mathematical operation used to combine the new spatial data with the geoprofile when operation = "continuous".
#' 
#' @export
#' @examples
#' # to come

geoMask <- function (probSurface, params, mask, scaleValue = 1, operation = "inside", maths = "multiply") {
  
  # check input formats
  stopifnot(class(mask) %in% c("SpatialPolygonsDataFrame", "RasterLayer","SpatialLinesDataFrame"))
  stopifnot(operation %in% c("inside", "outside", "near", "far", "continuous"))
  stopifnot(maths %in% c("multiply", "divide", "add", "subtract", "continuous"))
  
  # convert mask to raster
  if (class(mask) == "RasterLayer") { 
    rf <- mask
  } else if (class(mask) == "SpatialPolygonsDataFrame") {
    tmp <- raster::raster(ncol = params$output$longitude_cells, nrow = params$output$latitude_cells)
    raster::extent(tmp) <- raster::extent(mask)
    rf <- raster::rasterize(mask, tmp)
  } else if (class(mask) == "SpatialLinesDataFrame") {
    tmp <- raster(ncol = params$output$longitude_cells, nrow = params$output$latitude_cells)
    extent(tmp) <- extent(mask)
    rf <- rasterize(mask, tmp)
  }
  
  # convert probSurface to raster
  raster_probSurface <- raster(probSurface, xmn = params$output$longitude_minMax[1], xmx = params$output$longitude_minMax[2], ymn = params$output$latitude_minMax[1], ymx = params$output$latitude_minMax[2], crs="+proj=longlat +datum=WGS84")
  
  # project mask onto same coordinate system as probSurface
  rf <- raster::projectRaster(rf, raster_probSurface, crs="+proj=longlat +datum=WGS84")
  
  # extract raster values to matrices
  rf_mat <- matrix(rf@data@values, ncol = rf@ncols, byrow = TRUE)
  rf_mat <- rf_mat[nrow(rf_mat):1,]
  p_mat <- matrix(raster_probSurface@data@values, ncol = raster_probSurface@ncols, byrow = TRUE)
  
  #### OPERATIONS
  
  # initialise scale matrix
  scale_mat <- NULL
  
  # keep cells inside mask, multiplied by scaleValue
  if (operation == "inside") {
    scale_mat <- scaleValue * ifelse(is.na(rf_mat), NA, 1)
    p_mat <- p_mat * scale_mat
  }
  
  # keep cells outside mask, multiplied by scaleValue
  if (operation == "outside") {
    scale_mat <- scaleValue * ifelse(is.na(rf_mat), 1, NA)
    p_mat <- p_mat * scale_mat
  }
  
  # perform operation at all cells
  if (operation == "continuous") {
    if (maths == "add") { 
      p_mat <- p_mat + rf_mat
    }
    if (maths == "subtract") {
      p_mat <- p_mat - rf_mat
    }
    if (maths == "multiply") {
      p_mat <- p_mat * rf_mat
    }
    if (maths == "divide") {
      p_mat <- p_mat / rf_mat
    }
  }
  
  # decay with distance from non-NA cells
  if (operation == "near") {
    d <- raster::distance(rf)
    d_mat <- matrix(d@data@values, ncol = d@ncols, byrow = TRUE)
    scale_mat <- 1/(d_mat^scaleValue)
    scale_mat[scale_mat == "Inf"] <- 1
    p_mat <- p_mat * scale_mat
  }
  
  # increase with distance from non-NA cells
  if (operation == "far") {
    d <- raster::distance(rf)
    d_mat <- matrix(d@data@values, ncol = d@ncols, byrow = TRUE)
    scale_mat <- d_mat^scaleValue
    p_mat <- p_mat * scale_mat
  }
  
  # return list
  ret <- list(prob = p_mat, scaleMatrix = scale_mat)
  return(ret)
}
