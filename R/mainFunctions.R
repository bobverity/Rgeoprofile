# TO DO
# - the function ringHS() uses a lot of packages and seems quite complicated! At the moment all we did was change the unput variable names, but not the names in the function, meaning it won't currently work. If you could have a go at overhauling this function with the minimum dependencies that would be great. For example, if we can get away with adding our own circles rather than using spatialDataFrame objects that would be ideal.
# - expand help text and examples where needed. Remember that you need to run the function document() to actually create the help files once you've updated the text here. (NB, document() is part of devtools).
# - in functions geoReportHitscores() and geoPlotLorenz(), add comments and check program flow. Give input arguments defaults where possible (e.g. NULL), and do some formatting checks on inputs in case of bad input (e.g. using stopifnot())
# - separate functions that print results and those that return results, or just get rid of printing from existing functions. Make sure returned values are e.g. data frames with correctly labelled columns (particularly hitscores function)

#------------------------------------------------
#' Draw from Dirichlet process mixture model
#'
#' Provides random draws from a 2D spatial Dirichlet process mixture model. Both sigma and tau are defined in units of kilometres, representing the average distance that an observations lies from a source, and the average distance that a source lies from the centre point respectively. In contrast, the location of the centre point and the locations of the final output crime sites are defined in units of degrees lat/long to facilitate spatial analysis.
#'
#' Output includes the lat/long locations of the points drawn from the DPM model, along with the underlying group allocation (i.e. which points belong to which sources) and the lat/long locations of the sources.
#'
#' @param n number of draws.
#' @param sigma standard deviation of dispersal distribution, in units of kilometres.
#' @param tau standard deviation of prior on source locations (i.e. average distances of sources from centre point), in units of kilometres.
#' @param priorMean_longitude location of prior mean on source locations in degrees longitude.
#' @param priorMean_latitude location of prior mean on source locations in degrees latitude.
#' @param alpha concentration parameter of Dirichlet process model. Large alpha implies many distinct sources, while small alpha implies only a few sources.
#'
#' @export
#' @examples
#' # produces clusters of points from sources centred on QMUL
#' rDPM(50, priorMean_longitude = -0.04217491, priorMean_latitude = 51.5235505, 
#' alpha=1, sigma=1, tau=3) 
#' # same, but increasing alpha to generate more clusters
#' #' rDPM(50, priorMean_longitude = -0.04217491, priorMean_latitude = 51.5235505, 
#' alpha=5, sigma=1, tau=3)

rDPM <- function(n, sigma=1, tau=10, priorMean_longitude=-0.1277, priorMean_latitude=51.5074, alpha=1) {
	
	# force n to be a scalar integer
	n <- floor(n[1])
	
	# draw grouping
	group <- rep(1,n)
	freqs <- 1
    if (n>1) {
        for (i in 2:n) {
            group[i] <- sample(length(freqs)+1,1,prob=c(freqs,alpha))
            if (group[i]>length(freqs))
                freqs <- c(freqs,0)
            freqs[group[i]] <- freqs[group[i]] + 1
        }
    }
	group <- sort(group)
	
	# draw source locations
	source <- rnorm_sphere(length(freqs), priorMean_latitude, priorMean_longitude, tau)
	
	# draw crime locations from sources
	crime <- rnorm_sphere(n, source$latitude[group], source$longitude[group], sigma)
	
	# return results
	return(list(longitude=crime$longitude, latitude=crime$latitude, group=group, source_lon=source$longitude, source_lat=source$latitude))
}

#------------------------------------------------
# Calculate bearing and great circle distance between an origin and one or more destination points
# (not exported)

latlon_to_bearing <- function(origin_lat, origin_lon, dest_lat, dest_lon) {
	
	# convert input arguments to radians	
	origin_lat <- origin_lat*2*pi/360
	dest_lat <- dest_lat*2*pi/360
	origin_lon <- origin_lon*2*pi/360
	dest_lon <- dest_lon*2*pi/360
	
	delta_lon <- dest_lon-origin_lon
	
	# calculate bearing and great circle distance
	bearing <- atan2(sin(delta_lon)*cos(dest_lat), cos(origin_lat)*sin(dest_lat)-sin(origin_lat)*cos(dest_lat)*cos(delta_lon))
	gc_angle <- acos(sin(origin_lat)*sin(dest_lat) + cos(origin_lat)*cos(dest_lat)*cos(delta_lon))
	
	# convert bearing from radians to degrees measured clockwise from due north, and convert gc_angle to great circle distance via radius of earth (km)
	bearing <- bearing*360/(2*pi)
	bearing <- (bearing+360)%%360
	earthRad <- 6371
	gc_dist <- earthRad*gc_angle
	
	return(list(bearing=bearing, gc_dist=gc_dist))
}

#------------------------------------------------
# Convert lat/lon coordinates to cartesian coordinates by first calculating great circle distance and bearing and then mapping these coordinates into cartesian space. This mapping is relative to the point {centre_lat, centre_lon}, which should be roughly at the midpoint of the observed data.
# (not exported)

latlon_to_cartesian <- function(centre_lat, centre_lon, data_lat, data_lon) {
	
	# calculate bearing and great circle distance of data relative to centre
	data_trans <- latlon_to_bearing(centre_lat, centre_lon, data_lat, data_lon)
	
	# use bearing and distance to calculate cartesian coordinates
	theta <- data_trans$bearing*2*pi/360
	d <- data_trans$gc_dist
	data_x <- d*sin(theta)
	data_y <- d*cos(theta)
	
	return(list(x=data_x, y=data_y))
}

#------------------------------------------------
# Calculate destination lat/lon given an origin, a bearing and a great circle distance of travel
# Note that bearing should be in degrees relative to due north, and gc_dist should be in units of kilometres
# (not exported)

bearing_to_latlon <- function(origin_lat, origin_lon, bearing, gc_dist) {
	
	# convert origin_lat, origin_lon and bearing from degrees to radians
	origin_lat <- origin_lat*2*pi/360
	origin_lon <- origin_lon*2*pi/360
	bearing <- bearing*2*pi/360
	
	# calculate new lat/lon using great circle distance
	earthRad <- 6371
	new_lat <- asin(sin(origin_lat)*cos(gc_dist/earthRad) + cos(origin_lat)*sin(gc_dist/earthRad)*cos(bearing))
	new_lon <- origin_lon + atan2(sin(bearing)*sin(gc_dist/earthRad)*cos(origin_lat), cos(gc_dist/earthRad)-sin(origin_lat)*sin(new_lat))
	
	# convert new_lat and new_lon from radians to degrees
	new_lat <- new_lat*360/(2*pi)
	new_lon <- new_lon*360/(2*pi)
	
	return(list(longitude=new_lon, latitude=new_lat))
}

#------------------------------------------------
# Convert cartesian coordinates to lat/lon by using angle and euclidian distance from origin to define a bearing and great-circle distance relative to some centre point.
# (not exported)

cartesian_to_latlon <- function(centre_lat, centre_lon, data_x, data_y) {
	
	# calculate angle and euclidian distance of all points relative to origin
	d <- sqrt(data_x^2+data_y^2)
	theta <- atan2(data_y,data_x)
	
	# convert theta to bearing relative to due north
	theta <- theta*360/(2*pi)
	theta <- (90-theta)%%360
	
	# use bearing and great circle distance to calculate lat/lon relative to an origin point
	data_trans <- bearing_to_latlon(centre_lat, centre_lon, theta, d)
	
	return(list(longitude=data_trans$longitude, latitude=data_trans$latitude))
}

#------------------------------------------------
# Draw from normal distribution converted to spherical coordinate system. Points are first drawn from an ordinary cartesian 2D normal distribution. The distances to points are then assumed to be great circle distances, and are combined with a random bearing from the point {centre_lat, centre_lon} to produce a final set of lat/lon points. Note that this is not a truly spherical normal distribution, as the domain of the distribution is not the sphere - rather it is a transformation from one coordinate system to another that is satisfactory when the curvature of the sphere is not severe.
# (not exported)

rnorm_sphere <- function(n, centre_lat, centre_lon, sigma) {
	x <- rnorm(n,sd=sigma)
	y <- rnorm(n,sd=sigma)
	output <- cartesian_to_latlon(centre_lat, centre_lon, x, y)
	return(output)
}

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
#' geoData(Cholera[,1],Cholera[,2])
#' 
#' # simulated data
#' sim <-rDPM(50, priorMean_longitude = -0.04217491, priorMean_latitude = 
#' 51.5235505, alpha=1, sigma=1, tau=3)
#' geoData(sim$longitude, sim $latitude)

geoData <- function(longitude=NULL, latitude=NULL) {
	
	# use example data if none read in
	if (is.null(longitude) & is.null(latitude)) {
		data(LondonExample_crimes)
  	    longitude <- LondonExample_crimes$longitude
		latitude <- LondonExample_crimes$latitude
	} else {
		if (is.null(longitude) | is.null(latitude)) {
			stop("Both longitude and latitude arguments must be used, or alternatively both arguments must be set to NULL to use default values")
		}
	}
	
	# combine and return
	data <- list(longitude=longitude, latitude=latitude)
	return(data)
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
#' geoDataSource(WaterPumps[,1], WaterPumps[,2])
#' 
#' # simulated data
#' sim <-rDPM(50, priorMean_longitude = -0.04217491, priorMean_latitude = 
#' 51.5235505, alpha=1, sigma=1, tau=3)
#' geoDataSource(sim$source_longitude, sim$source_latitude)

geoDataSource <- function(source_longitude=NULL, source_latitude=NULL) {
    
  # generate dummy data if none read in
  if (is.null(source_longitude) & is.null(source_latitude)) {
		data(LondonExample_sources)
  	    source_longitude <- LondonExample_sources$longitude
		source_latitude <- LondonExample_sources$latitude
  } else {
  	if (is.null(source_longitude) | is.null(source_latitude))
  		stop("Both source_longitude and source_latitude arguments must be used, or alternatively both arguments must be set to NULL to use default values")
  }
  
  # combine and return
  source_data <- list(source_longitude=source_longitude, source_latitude=source_latitude)
  return(source_data)
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
#' d <- geoData(Cholera[,1], Cholera[,2])
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
        if (is.null(priorMean_longitude))
        	priorMean_longitude <- sum(range(data$longitude))/2
        if (is.null(priorMean_latitude))
        	priorMean_latitude <- sum(range(data$latitude))/2        
        
        # convert data to bearing and great circle distance, and extract maximum great circle distance to any point
        data_trans <- latlon_to_bearing(priorMean_latitude, priorMean_longitude, data$latitude, data$longitude)
		dist_max <- max(data_trans$gc_dist)
		
       	# use maximum distance as default value of tau
       	if (is.null(tau))
       		tau <- dist_max
		
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
	
		if (is.null(priorMean_longitude))
        	priorMean_longitude <- -0.1277
		if (is.null(priorMean_latitude))
        	priorMean_latitude <- 51.5074
		if (is.null(longitude_minMax))
			longitude_minMax <- priorMean_longitude + c(-0.01,0.01)
		if (is.null(latitude_minMax))
			latitude_minMax <- priorMean_latitude + c(-0.01,0.01)
       	if (is.null(tau))
       		tau <- 10
    }
	
    # initialise shape and rate parameters for prior on sigma^2
    alpha <- NULL
    beta <- NULL
    
    # if sigma_var has been specified
    if (!is.null(sigma_var)) {
    	
    	# check that sigma_mean has also been specified
    	if (!is.null(sigma_mean)) {
    		
    		# if using fixed sigma model then no need to calculate alpha and beta. Otherwise use values of sigma_mean and sigma_var to search for the unique alpha and beta that define the distribution
	    	if (sigma_var==0) {
	    		cat('Using fixed sigma model')
	    	} else {
	            cat('Using sigma_mean and sigma_var to define prior on sigma')
	    		ab <- get_alpha_beta(sigma_mean, sigma_var)
		    	alpha <- ab$alpha
		    	beta <- ab$beta
	    	}
    	}
    	    	
    # if sigma_var has not been specified
    } else {
    	
    	# if sigma_mean has been specified but sigma_var has not then use sigma_mean along with sigma_squared_shape to calculate beta
    	if (!is.null(sigma_mean)) {
	        if (is.null(sigma_squared_shape)) {
	        	stop("Current prior parameters on sigma do not fully specify the distribution. Must specify either 1) a prior mean and variance on sigma, 2) a prior mean on sigma and a prior shape on sigma^2, 3) a prior shape and prior rate on sigma^2.")
	        } else {
	            cat('Using sigma_mean and sigma_squared_shape to define prior on sigma')
	        	alpha <- sigma_squared_shape
	        	if (alpha<=1) {
	        		stop('sigma_squared_shape must be >1')
	        	}
	            beta <- exp(2*log(sigma_mean) + 2*lgamma(alpha) - 2*lgamma(alpha-0.5))
	            epsilon <- sqrt(beta)*gamma(alpha-0.5)/gamma(alpha)
				sigma_var <- beta/(alpha-1)-epsilon^2
	        }
	    
	    # if neither sigma_mean nor sigma_var have been specified then use sigma_squared_shape and sigma_squared_rate to define distribution
	    } else {
	    	cat('Using sigma_squared_shape and sigma_squared_rate to define prior on sigma')
	    	if (is.null(sigma_squared_shape) | is.null(sigma_squared_rate)) {
				stop("Current prior parameters on sigma do not fully specify the distribution. Must specify either 1) a prior mean and variance on sigma, 2) a prior mean on sigma and a prior shape on sigma^2, 3) a prior shape and prior rate on sigma^2.")
	    	}
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
    params <- list(model=model, MCMC=MCMC, output=output)
    return(params)
}

#------------------------------------------------
# Get alpha and beta parameters of inverse-gamma prior on sigma^2 from expectation and variance.
# (not exported)

get_alpha_beta <- function(sigma_mean,sigma_var) {
  
  # define a function that has minimum at correct value of alpha
  f_alpha <- function(alpha) {
    (sqrt((sigma_var+sigma_mean^2)*(alpha-1))*exp(lgamma(alpha-0.5)-lgamma(alpha))-sigma_mean)^2
  }
  
  # search for alpha
  alpha <- optim(2,f_alpha,method='Brent',lower=1,upper=1e3)$par
  
  # solve for beta
  beta <- (sigma_var+sigma_mean^2)*(alpha-1)
  
  # check that chosen alpha is not at limit of range
  if (alpha>(1e3-1))
      stop('unable to define prior on sigma for chosen values of sigma_mean and sigma_var. Try increasing the value of sigma_var, or alternatively setting sigma_var=0 (i.e. using fixed-sigma model)')
      
  output <- list(alpha=alpha, beta=beta)
  return(output)
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
#' d <- geoData(Cholera[,1],Cholera[,2])
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
  if (!("longitude"%in%names(data)))
    stop("data must contain element 'longitude'")
  if (!("latitude"%in%names(data)))
    stop("data must contain element 'latitude'")
  
  # check that data values are correct format and range
  if (!is.numeric(data$longitude) | !all(is.finite(data$longitude)))
    stop("data$longitude values must be numeric and finite")
  if (!is.numeric(data$latitude) | !all(is.finite(data$latitude)))
    stop("data$latitude values must be numeric and finite")
  
  # check same number of observations in logitude and latitude, and n>1
  if (length(data$longitude)!=length(data$latitude))
    stop("data$longitude and data$latitude must have the same length")
  if (length(data$longitude)<=1)
    stop("data must contain at least two observations")
  
  # if passed all checks
  if (!silent)
    cat("data file passed all checks\n")
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
#' d <- geoData(Cholera[,1],Cholera[,2])
#' p <- geoParams(data = d, sigma_mean = 1.0, sigma_squared_shape = 2)
#' geoParamsCheck(p)
#' 
#' # simulated data
#' sim <-rDPM(50, priorMean_longitude = -0.04217491, priorMean_latitude = 
#' 51.5235505, alpha=1, sigma=1, tau=3)
#' d <- geoData(Cholera[,1], Cholera[,2])
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
  if (!silent)
    cat("params file passed all checks\n")
}

#------------------------------------------------
#' Plot prior and posterior distributions of sigma.
#'
#' Plot prior distribution of sigma as defined by current parameter values. Can optionally overlay a kernel density plot of posterior draws of sigma.
#'
#' @param params a list of parameters as defined by geoParams().
#' @param mcmc stored output obtained by running geoMCMC(). Leave as NULL to plot prior only.
#' @param plotMax maximum x-axis range to plot. Leave as NULL to use default settings.
#'
#' @export
#' @examples
#' # John Snow cholera data
#' data(Cholera)
#' d <- geoData(Cholera[,1],Cholera[,2])
#' p <- geoParams(data = d, sigma_mean = 0.2, sigma_squared_shape = 2)
#' m <- geoMCMC(data = d, params = p)
#' geoPlotSigma(params = p, mcmc = m)
#' 
#' # simulated data
#' sim <-rDPM(50, priorMean_longitude = -0.04217491, priorMean_latitude = 
#' 51.5235505, alpha=1, sigma=1, tau=3)
#' d <- geoData(sim$longitude, sim$latitude)
#' p <- geoParams(data = d, sigma_mean = 1.0, sigma_squared_shape = 2)
#' m <- geoMCMC(data = d, params = p)
#' geoPlotSigma(params = p, mcmc = m)

geoPlotSigma <- function(params, mcmc=NULL, plotMax=NULL) {
  
  # check params
  geoParamsCheck(params, silent=TRUE)
  
  # check that plotMax is sensible
  if (!is.null(plotMax)) {
    if (!is.numeric(plotMax) | !is.finite(plotMax))
      stop('plotMax must be numeric and finite')
    if (plotMax<=0)
      stop('plotMax must be greater than zero')
  }
  
  # extract sigma parameters
  sigma_mean <- params$model$sigma_mean
  sigma_var <- params$model$sigma_var
  alpha <- params$model$sigma_squared_shape
  beta <- params$model$sigma_squared_rate
  
  # stop if using fixed sigma model
  if (sigma_var==0)
    stop('can only produce this plot under variable-sigma model (i.e. sigma_var>0)')
  
  # extract sigma draws from mcmc object
  sigma_draws <- mcmc$sigma
  
  # default plotMax based on extent of prior distribution AND the extent of posterior draws if available
  if (is.null(plotMax)) {
    plotMax <- sigma_mean+3*sqrt(sigma_var)
    if (!is.null(sigma_draws)) {
      plotMax <- max(plotMax, 2*max(sigma_draws,na.rm=TRUE))
    }
  }
  
  # produce prior distribution
  sigma_vec <- seq(0,plotMax,l=501)
  sigma_prior <- dRIG(sigma_vec,alpha,beta)
  
  # plot prior and overlay density of posterior draws if used
  if (is.null(sigma_draws)) {
  	
    plot(sigma_vec, sigma_prior, type='l', xlab='sigma (km)', ylab='probability density', main='')
    legend(x='topright', legend='prior', lty=1)
    
  } else {
  	
    sigma_posterior <- density(sigma_draws,from=0,to=plotMax)
    y_max <- max(sigma_posterior$y,na.rm=TRUE)
    
    plot(sigma_vec, sigma_prior, type='l', ylim=c(0,y_max), lty=2, xlab='sigma (km)', ylab='probability density', main='')
    lines(sigma_posterior)
    legend(x='topright', legend=c('prior','posterior'), lty=c(2,1))
  }
  
}

#------------------------------------------------
# Square-root-inverse-gamma distribution
# If an inverse gamma distribution has shape alpha and rate beta, and hence mean beta/(alpha-1) and variance beta^2/((alpha-1)^2*(alpha-2)), then the square root of this random variable has mean epsilon=sqrt(beta)*gamma(alpha-0.5)/gamma(alpha) and variance v=beta/(alpha-1)-epsilon^2. The variance can also be written purely in terms of alpha and epsilon as follows: v=epsilon^2*(gamma(alpha-1)*gamma(alpha)/gamma(alpha-0.5)^2 - 1).
# (not exported)

dRIG <- function(x,alpha,beta,log=FALSE) {
  output <- log(2)+alpha*log(beta)-lgamma(alpha)-(2*alpha+1)*log(x)-beta/x^2
  if (!log)
    output <- exp(output)
  return(output)
}

#------------------------------------------------
# Scaled Student's t distribution. Used in kernel density smoothing.
# (not exported)

dts <- function(x,df,scale=1,log=FALSE) {
  output <- lgamma((df+1)/2)-lgamma(df/2)-0.5*log(pi*df*scale^2) - ((df+1)/2)*log(1 + x^2/(df*scale^2))
  if (!log)
    output <- exp(output)
  return(output)
}

#------------------------------------------------
# Bin values in two dimensions
# (not exported)

bin2D <- function(x, y, x_breaks, y_breaks) {
    
    # get number of breaks in each dimension
    nx <- length(x_breaks)
    ny <- length(y_breaks)
    
    # create table of binned values
    tab1 <- table(findInterval(x, x_breaks), findInterval(y, y_breaks))
    
    # convert to dataframe and force numeric
    df1 <- as.data.frame(tab1, stringsAsFactors=FALSE)
    names(df1) <- c("x", "y", "count")
    df1$x <- as.numeric(df1$x)
    df1$y <- as.numeric(df1$y)
    
    # subset to within breaks range
    df2 <- subset(df1, x>0 & x<nx & y>0 & y<ny)
    
    # fill in matrix
    mat1 <- matrix(0,ny-1,nx-1)
    mat1[cbind(df2$y, df2$x)] <- df2$count
    
    # calculate cell midpoints
    x_mids <- (x_breaks[-1]+x_breaks[-nx])/2
    y_mids <- (y_breaks[-1]+y_breaks[-ny])/2
    
    # return output as list
    output <- list(x_mids=x_mids, y_mids=y_mids, z=mat1)
    return(output)
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
#' d <- geoData(Cholera[,1],Cholera[,2])
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
    
    # bin mu draws in two dimensions and check that at least one posterior draw in chosen region
    surface_raw <- bin2D(mu_draws$longitude, mu_draws$latitude, breaks_lon, breaks_lat)$z
    if (all(surface_raw==0))
    stop('chosen lat/long window contains no posterior draws')
    
    # temporarily add guard rail to surface to avoid Fourier series bleeding round edges
    railSize_x <- cells_lon
    railSize_y <- cells_lat
    railMat_x <- matrix(0,cells_lat,railSize_x)
    railMat_y <- matrix(0,railSize_y,cells_lon+2*railSize_x)
    
    surface_normalised <- surface_raw/sum(surface_raw)
    surface_normalised <- cbind(railMat_x, surface_normalised, railMat_x)
    surface_normalised <- rbind(railMat_y, surface_normalised, railMat_y)
    
    # calculate Fourier transform of posterior surface
    f1 = fftw2d(surface_normalised)
    
    # calculate x and y size of one cell in cartesian space. Because of transformation this size will technically be different for each cell, but use centre of prior to get a middling values
    cellsize_prior <-latlon_to_cartesian(params$model$priorMean_latitude, params$model$priorMean_longitude, params$model$priorMean_latitude+cellSize_lat, params$model$priorMean_longitude+cellSize_lon)
    cellsize_x <- cellsize_prior$x
    cellsize_y <- cellsize_prior$y
    
    # produce surface over which kernel will be calculated. This surface wraps around in both x and y (i.e. the kernel is actually defined over a torus).
    kernel_x <- cellsize_x * c(0:floor(ncol(surface_normalised)/2), floor((ncol(surface_normalised)-1)/2):1)
    kernel_y <- cellsize_y * c(0:floor(nrow(surface_normalised)/2), floor((nrow(surface_normalised)-1)/2):1)
    kernel_x_mat <- outer(rep(1,length(kernel_y)), kernel_x)
    kernel_y_mat <- outer(kernel_y, rep(1,length(kernel_x)))
    kernel_s_mat <- sqrt(kernel_x_mat^2+kernel_y_mat^2)
    
    # set lambda (bandwidth) range to be explored
    if (is.null(lambda)) {
	    lambda_step <- min(cellsize_x, cellsize_y)/5
	    lambda_vec <- lambda_step*(1:100)
	} else {
		lambda_vec <- lambda
	}
    
    # loop through range of values of lambda
    cat('Smoothing posterior surface')
    flush.console()
    logLike <- -Inf
    for (i in 1:length(lambda_vec)) {
        
        if (i>1) {
        	cat(".")
        	flush.console()
        }
        
        # calculate Fourier transform of kernel
        lambda_this <- lambda_vec[i]
        kernel <- dts(kernel_s_mat,df=3,scale=lambda_this)
        f2 = fftw2d(kernel)
        
        # combine Fourier transformed surfaces and take inverse. f4 will ultimately become the main surface of interest.
        f3 = f1*f2
        f4 = Re(fftw2d(f3,inverse=T))/length(surface_normalised)
        
        # subtract from f4 the probability density of each point measured from itself. In other words, move towards a leave-one-out kernel density method
        f5 <- f4 - surface_normalised*dts(0,df=3,scale=lambda_this)
        f5[f5<0] <- 0
        f5 <- f5/sum(f4)
        
        # calculate leave-one-out log-likelihood at each point on surface
        f6 <- surface_normalised*log(f5)
        
        # break if total log-likelihood is at a local maximum
        if (sum(f6,na.rm=T)<logLike) {
	        break()
	    }
        logLike <- sum(f6,na.rm=T)
        
    }
    
    # report chosen value of lambda
    if (is.null(lambda)) {
	    cat(paste('\nmaximum likelihood lambda = ',round(lambda_this,3),sep=''))
	}
    
    # remove guard rail
    f4 <- f4[,(railSize_x+1):(ncol(f4)-railSize_x)]
    f4 <- f4[(railSize_y+1):(nrow(f4)-railSize_y),]
    
    # calculate coordinates of lat/lon matrix in original cartesian coordinates
    cart <-latlon_to_cartesian(params$model$priorMean_latitude, params$model$priorMean_longitude, mids_lat_mat, mids_lon_mat)
    
    # produce prior matrix. Note that each cell of this matrix contains the probability density at that point multiplied by the size of that cell, meaning the total sum of the matrix from -infinity to +infinity would equal 1. However, as the matrix is limited to the region specified by the limits, in reality this matrix will sum to some value less than 1.
    priorMat <- dnorm(cart$x,sd=params$model$tau)*dnorm(cart$y,sd=params$model$tau)*(cellSize_lon*cellSize_lat)
    
    # combine prior surface with stored posterior surface (the prior never fully goes away under a DPM model)
    n <- length(data$longitude)
    alpha <- rawOutput$alpha
    posteriorMat <-  f4 + priorMat*mean(alpha/(alpha+n))
    
    # produce geoprofile
    gp <- geoProfile(posteriorMat)
    
    # calculate posterior allocation
    allocation <- matrix(unlist(rawOutput$allocation),n,byrow=T)
    allocation <- data.frame(allocation/params$MCMC$samples)
    names(allocation) <- paste("group",1:ncol(allocation),sep="")
    
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
    
    return(output)
}

#------------------------------------------------
#' Calculate geoprofile from surface
#'
#' Converts surface to rank order geoprofile.
#'
#' @param surface matrix to convert to geoprofile
#'
#' @export
#' @examples
#' # John Snow cholera data
#' data(Cholera)
#' d <- geoData(Cholera[,1],Cholera[,2])
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
    if (!all(is.finite(surface)) | !all(is.numeric(surface)))
        stop("values in surface must be finite and numeric")
    if (!is.matrix(surface))
        stop("surface must be a matrix")
    
    # create geoprofile from surface
    surface[order(surface,decreasing=TRUE)] <- 1:length(surface)
    
    return(surface)
}

#------------------------------------------------
#' Plot posterior allocation
#'
#' Produces plot of posterior allocation from output of MCMC.
#'
#' @param mcmc stored output obtained by running geoMCMC().
#' @param colors vector of colours for each allocation. If NULL then use default colour scheme.
#' @param barBorderCol colour of borders around each bar. Set as NA to omit this border (useful when there are a large number of observations).
#' @param barBorderWidth line width of borders around each bar.
#' @param mainBorderCol colour of border around plot.
#' @param mainBorderWidth line width of border around plot.
#' @param yTicks_on whether to include ticks on the y-axis.
#' @param yTicks vector of y-axis tick positions.
#' @param xTicks_on whether to include ticks on the x-axis.
#' @param xTicks_size size of ticks on the x-axis.
#' @param xlab x-axis label.
#' @param ylab x-axis label.
#' @param mainTitle main title over plot.
#' @param names individual names of each observation, written horizontally below each bar.
#' @param names_size size of names under each bar.
#' @param orderBy whether to order segments within each bar by "group" or by "probability". If ordered by group, all segments of a particular group are laid down before moving to the next group. If ordered by probability the segments within each bar are ordered from large to small.
#'
#' @export
#' @examples
#' # London example data
#' d <- geoData()
#' p <- geoParams(data = d, sigma_mean = 1.0, sigma_squared_shape = 2)
#' m <- geoMCMC(data = d, params = p)
#' geoPlotAllocation(m)
#'
#' # John Snow cholera data
#' data(Cholera)
#' d <- geoData(Cholera[,1],Cholera[,2])
#' p <- geoParams(data = d, sigma_mean = 1.0, sigma_squared_shape = 2)
#' m <- geoMCMC(data = d, params = p, lambda=0.05)
#' geoPlotAllocation(m, barBorderCol=NA)	# (should allocate all to a single source!)
#' 
#' # simulated data
#' sim <-rDPM(50, priorMean_longitude = -0.04217491, priorMean_latitude = 
#' 51.5235505, alpha=1, sigma=1, tau=3)
#' d <- geoData(sim$longitude, sim $latitude)
#' p <- geoParams(data = d, sigma_mean = 1.0, sigma_squared_shape = 2)
#' m <- geoMCMC(data = d, params = p)
#' geoPlotAllocation(m)

geoPlotAllocation <- function(mcmc, colours="default", barBorderCol="white", barBorderWidth=0.25, mainBorderCol="black", mainBorderWidth=2, yTicks_on=TRUE, yTicks=seq(0,1,0.2), xTicks_on=FALSE, xTicks_size=1, xlab="", ylab="posterior allocation", mainTitle="", names=NA, names_size=1, orderBy="group") {
    
    # check that orderBy is either 'group' or 'probability'
    if (!(orderBy%in%c("group","probability")))
        stop("orderBy must equal 'group' or 'probability'")
    
    # get allocation from mcmc object
    allocation <- mcmc$allocation
    
    # check that allocation is a data frame
    if (!is.data.frame(allocation))
        stop("allocation must be a data frame, with observations in rows and groups in columns")
    
    n <- nrow(allocation)
    k <- ncol(allocation)
    
    # replace colours if default
    if (colours=="default") {
        rawCols <- brewer.pal(n=11,name="RdYlBu")
        myPal <- colorRampPalette(rawCols)
        colours <- myPal(k)
        colours <- colours[c(2*(1:ceiling(k/2))-1,2*(1:floor(k/2)))]
    }
    
    # if ordered by group
    if (orderBy=="group") {
        
        plot(0, type='n', xlim=c(0,n), ylim=c(0,1), xlab=xlab, ylab=ylab, xaxs="i", yaxs="i",axes=FALSE, main=mainTitle)
        barplot(t(allocation), col=colours, space=0, border=NA, axes=FALSE, add=TRUE)
        segments(0:n,c(0,0),0:n,c(1,1), col=barBorderCol, lwd=barBorderWidth)
        box(col=mainBorderCol, lwd=mainBorderWidth)
        axis(2, tick=yTicks_on, labels=yTicks_on, at=yTicks)
        axis(1, at=1:n-0.5, tick=xTicks_on, lwd.ticks=xTicks_size, labels=names, las=2, cex.axis=names_size)
    }
    
    # if ordered by probability
    if (orderBy=="probability") {
        
        plot(0, type='n', xlim=c(0,n), ylim=c(0,1), xlab=xlab, ylab=ylab, xaxs="i", yaxs="i",axes=FALSE, main=mainTitle)
        tM <- t(allocation)
        for (i in 1:n) {
            temp <- tM
            temp[,-i] <- NA
            temp_order <- order(temp[,i],decreasing=TRUE)
            barplot(temp[temp_order,], col=colours[temp_order], space=0, border=NA, axes=FALSE, add=TRUE)
        }
        segments(0:n,c(0,0),0:n,c(1,1), col=barBorderCol, lwd=barBorderWidth)
        box(col=mainBorderCol, lwd=mainBorderWidth)
        axis(2, tick=yTicks_on, labels=yTicks_on, at=yTicks)
        axis(1, at=1:n-0.5, tick=xTicks_on, lwd.ticks=xTicks_size, labels=names, las=2, cex.axis=names_size)
    }
    
}

#------------------------------------------------
# get optimal zoom level given x and y values
# (not exported)

getZoom <- function(x,y) {
	
	# calculate midpoint of range in x and y
	xmid <- min(x)+diff(range(x))/2
	ymid <- min(y)+diff(range(y))/2
	
	# calculate delta (half of longitude range) for a range of zoom levels
	z <- 20:2
	delta <- 445/(2^z)
	
	# calculate left and right longitude limits at all zoom levels
	long_angle_left <- xmid-delta
	long_angle_right <- xmid+delta
	
	# calculate top and bottom latitude limits at all zoom levels
	lat_angle <- ymid/360*2*pi
	projection_mid <- log(tan(pi/4+lat_angle/2))
	projection_top <- projection_mid + delta/360*2*pi
	projection_bot <- projection_mid - delta/360*2*pi
	lat_angle_top <- (2*atan(exp(projection_top))-pi/2)*360/(2*pi)
	lat_angle_bot <- (2*atan(exp(projection_bot))-pi/2)*360/(2*pi)
	
	# find the most zoomed-in level that captures all points
	zoomTest <- (min(x)>long_angle_left) & (max(x)<long_angle_right) & (min(y)>lat_angle_bot) & (max(y)<lat_angle_top)
    if (!any(zoomTest))
        stop("values are outside of plotting range of earth")
	bestZoom <- z[which(zoomTest)[1]]
    
	return(bestZoom)
}

#------------------------------------------------
#' Plot a map and overlay data and/or geoprofile
#'
#' Plots geoprofile on map, with various customisable options.
#'
#' @param params parameters list in the format defined by geoParams().
#' @param data data object in the format defined by geoData().
#' @param source potential sources object in the format defined by geoDataSource().
#' @param surface a surface to overlay onto the map, typically a geoprofile obtained from the output of geoMCMC().
#' @param zoom zoom level of map. If NULL then choose optimal zoom from params.
#' @param mapSource which online source to use when downloading the map. Options include Google Maps ("google"), OpenStreetMap ("osm"), Stamen Maps ("stamen") and CloudMade maps ("cloudmade").
#' @param mapType the specific type of map to plot. Options available are "terrain", "terrain-background", "satellite", "roadmap" and "hybrid" (google maps), "terrain", "watercolor" and "toner" (stamen maps) or a positive integer for cloudmade maps (see ?get_cloudmademap from the package ggmap for details).
#' @param opacity value between 0 and 1 givin the opacity of surface colours.
#' @param plotContours whether or not to add contours to the surface plot.
#' @param breakPercent vector of values between 0 and 100 describing where in the surface contours appear.
#' @param contourCols list of two or more colours from which to derive the contour colours.
#' @param crimeCex relative size of symbols showing crimes.
#' @param crimeCol colour of crime symbols.
#' @param crimeBorderCol border colour of crime symbols.
#' @param crimeBorderWidth width of border of crime symbols.
#' @param sourceCex relative size of symbols showing suspect sites.
#' @param sourceCol colour of suspect sites symbols.
#' @param gpLegend whether or not to add legend to plot.
#'
#' @export
#' @examples
#' # London example data
#' d <- geoData()
#' s <- geoDataSource()
#' p <- geoParams(data = d, sigma_mean = 1.0, sigma_squared_shape = 2)
#' m <- geoMCMC(data = d, params = p)
#' # produce simple map
#' geoPlotMap(params = p, data = d, source = s, surface = m$geoProfile, breakPercent = seq(0, 50, 5), mapType = "hybrid",
#' crimeCol = "black", crimeCex = 2, sourceCol = "red", sourceCex = 2)
#'
#' # John Snow cholera data
#' data(Cholera)
#' d <- geoData(Cholera[,1],Cholera[,2])
#' data(WaterPumps)
#' s <- geoDataSource(WaterPumps[,1], WaterPumps[,2])
#' p <- geoParams(data = d, sigma_mean = 1.0, sigma_squared_shape = 2)
#' m <- geoMCMC(data = d, params = p, lambda=0.05)
#' # produce simple map
#' geoPlotMap(params = p, data = d, source = s, surface = m$geoProfile, breakPercent = seq(0, 50, 5), mapType = "hybrid",
#' crimeCol = "black", crimeCex = 2, sourceCol = "red", sourceCex = 2)
#' 
#' # simulated data
#' sim <-rDPM(50, priorMean_longitude = -0.04217491, priorMean_latitude = 
#' 51.5235505, alpha=1, sigma=1, tau=3)
#' d <- geoData(sim$longitude, sim $latitude)
#' s <- geoDataSource(sim$source_lon, sim$source_lat)
#' p <- geoParams(data = d, sigma_mean = 1.0, sigma_squared_shape = 2)
#' m <- geoMCMC(data = d, params = p)
#' # change colour palette, map type, opacity and range of geoprofile and omit legend
#' geoPlotMap(params = p, data = d, source = s, surface = m$geoProfile, breakPercent = seq(0, 30, 5), mapType = "terrain", 
#' contourCols=c("blue","white"),crimeCol="black", crimeBorderCol="white",crimeCex=2,
#' sourceCol = "red", sourceCex = 2, opacity = 0.7, gpLegend = FALSE)

geoPlotMap <- function(params, data=NULL, source=NULL, surface=NULL, zoom=NULL, mapSource="google", mapType="hybrid", opacity=0.6, plotContours=TRUE, breakPercent=seq(0,100,l=11), contourCols= c("red","orange","yellow","white"), crimeCex=1.5, crimeCol='red', crimeBorderCol='white', crimeBorderWidth=0.5, sourceCex=1.5, sourceCol='blue', gpLegend=TRUE) {
    
    # check that inputs make sense
    geoParamsCheck(params)
    if (!is.null(data))
	    geoDataCheck(data)
    
    # if zoom=="auto" then set zoom level based on params
    if (is.null(zoom))
        zoom <- getZoom(params$output$longitude_minMax, params$output$latitude_minMax)
    
    # make zoom level appropriate to map source
    if (mapSource=="stamen")
    	zoom <- min(zoom,18)
    
    # make map
    rawMap <- get_map(location=c(mean(params$output$longitude_minMax), mean(params$output$latitude_minMax)), zoom=zoom, source=mapSource, maptype=mapType)
    
    # add limits
    myMap <- ggmap(rawMap) + coord_cartesian(xlim=params$output$longitude_minMax, ylim=params$output$latitude_minMax)
    
    # overlay geoprofile
    if (!is.null(surface)) {
    	
    	# create colour palette
    	geoCols <- colorRampPalette(contourCols)
    	nbcol=length(breakPercent)-1
    	
    	# extract plotting ranges and determine midpoints of cells
    	longitude_minMax  <- params$output$longitude_minMax
    	latitude_minMax  <- params$output$latitude_minMax
    	longitude_cells  <- params$output$longitude_cells
    	latitude_cells  <- params$output$latitude_cells
    	longitude_cellSize <- diff(longitude_minMax)/longitude_cells
    	latitude_cellSize <- diff(latitude_minMax)/latitude_cells
    	longitude_midpoints <- longitude_minMax[1] - longitude_cellSize/2 + (1:longitude_cells)* longitude_cellSize
    latitude_midpoints <- latitude_minMax[1] - latitude_cellSize/2 + (1:latitude_cells)* latitude_cellSize
    	
    	# create data frame of x,y,z values and labels for contour level
		df <- expand.grid(x=longitude_midpoints, y=latitude_midpoints)
		df$z <- as.vector(t(surface))
		labs <- paste(round(breakPercent,1)[-length(breakPercent)],"-",round(breakPercent,1)[-1],"%",sep='')
		df$cut <- cut(df$z, breakPercent/100*length(surface), labels=labs)
		
		# remove all entries outside of breakPercent range
		df_noNA <- df[!is.na(df$cut),]
		
		# add surface and hitscore legend
		myMap <- myMap + geom_tile(aes(x=x,y=y,fill=cut), alpha=opacity, data=df_noNA)
		myMap <- myMap + scale_fill_manual(name="Hitscore\npercentage", values=rev(geoCols(nbcol)))
		if(gpLegend==FALSE) {myMap <- myMap + theme(legend.position="none")}

		# add contours
		if (plotContours) {
			myMap <- myMap + stat_contour(aes(x=x,y=y,z=z), colour="grey50", breaks=breakPercent/100*length(surface), size=0.3, alpha=opacity, data=df)
		}
	}

    # overlay data points
    if (!is.null(data)) {
    	df_data <- data.frame(longitude=data$longitude, latitude=data$latitude)
		myMap <- myMap + geom_point(aes(x=longitude, y=latitude), data=df_data, pch=21, stroke=crimeBorderWidth, cex=crimeCex, fill=crimeCol, col=crimeBorderCol)
    }
    
    # overlay source points
    if (!is.null(source)) {
    	df_source <- data.frame(longitude=source$source_longitude, latitude=source$source_latitude)
		myMap <- myMap + geom_point(aes(x=longitude, y=latitude), data=df_source, pch=15, cex=sourceCex, col=sourceCol)
    }
    
    # plot map
    myMap
}

#------------------------------------------------
#' Calculate hitscores
#'
#' Calculate hitscores of the potential sources for a given surface (usually the geoprofile).
#'
#' @param mcmc stored output obtained by running geoMCMC().
#' @param source longitude and latitude of one or more source locations in the format defined by geoDataSource().
#'
#' @export
#' @examples
#' # John Snow cholera data
#' data(Cholera)
#' d <- geoData(Cholera[,1],Cholera[,2])
#' data(WaterPumps)
#' s <- geoDataSource(WaterPumps[,1], WaterPumps[,2])
#' p <- geoParams(data = d, sigma_mean = 1.0, sigma_squared_shape = 2, samples = 20000, 
#' chains = 10, burnin = 1000, priorMean_longitude = mean(d$longitude), 
#' priorMean_latitude = mean(d$latitude), guardRail = 0.1)
#' m <- geoMCMC(data = d, params = p)
#' gp <- m$geoProfile
#' geoPlotMap(params = p, data = d, source = s, breakPercent = seq(0, 50, 5), mapType = "hybrid",
#' contourCols = c("red", "orange", "yellow", "white"), crimeCol = "black", crimeBorderCol = "white", 
#' crimeCex = 2, sourceCol = "red", sourceCex = 2, surface = gp)
#' hs <- geoReportHitscores(params=p,source_data=s,surface=m$surface)
#' 
#' # simulated data
#' sim <-rDPM(50, priorMean_longitude = -0.04217491, priorMean_latitude = 
#' 51.5235505, alpha=1, sigma=1, tau=3)
#' d <- geoData(sim$longitude, sim $latitude)
#' s <- geoDataSource(sim$source_lon, sim$source_lat)
#' p <- geoParams(data = d, sigma_mean = 1.0, sigma_squared_shape = 2, samples = 20000, 
#' chains = 10, burnin = 1000, priorMean_longitude = mean(d$longitude), 
#' priorMean_latitude = mean(d$latitude), guardRail = 0.1)
#' gp <- m$geoProfile
#' geoPlotMap(params = p, data = d, source = s,breakPercent = seq(0, 30, 5), mapType = "terrain", 
#' contourCols=c("blue","white"),crimeCol="black", crimeBorderCol="white",crimeCex=2,
#' sourceCol = "red", sourceCex = 2, surface = gp, transparency = 0.7)
#' hs <- geoReportHitscores(params=p,source_data=s,surface=m$surface)
#' 
geoReportHitscores <- function(params, source_data, surface) {
sources <- cbind(source_data$source_longitude,source_data$source_latitude)
	ordermat = matrix(0,params$output$latitude_cells,params$output$longitude_cells)
	profile_order = order(surface)
	for (i in 1:(params$output$latitude_cells * params$output$longitude_cells)) {
		ordermat[profile_order[i]] = i
		}
	hitscoremat <<- 1-ordermat/(params$output$latitude_cells * params$output$longitude_cells)
	hitscoremat2 <- hitscoremat[nrow(hitscoremat):1,]

	xvec=seq(params$output$longitude_minMax[1],params$output$longitude_minMax[2],length=params$output$longitude_cells)
	yvec=seq(params$output$latitude_minMax[1],params$output$latitude_minMax[2],length=params$output$latitude_cells)
xdiff = abs(outer(rep(1,nrow(sources)),xvec)-outer(sources[,1],rep(1,params$output$longitude_cells)))
ydiff = abs(outer(rep(1,nrow(sources)),yvec)-outer(sources[,2],rep(1,params$output$latitude_cells)))

msourcex = mapply(which.min,x=split(xdiff,row(xdiff)))
msourcey = params$output$longitude_cells-(mapply(which.min,x=split(ydiff,row(ydiff))))+1

	if (nrow(sources)>1) {
		hitscores = diag(hitscoremat2[msourcey,msourcex])
	} else {
		hitscores = hitscoremat2[msourcey,msourcex]}
hit_output <<- cbind(sources,hitscores)
	colnames(hit_output) <- c("lon","lat","hs")
	return(hit_output)
}
#------------------------------------------------
#' Produce Lorenz Plot
#'
#' Produces a Lorenz plot showing the proportion of suspect sites or cimes identified as a function of area and calculates
#' the corresponding Gini coefficient using trapezoid rule.
#' Also allows an optional vector called crimeNumbers with numbers of crimes per suspect site; tthe length of this vector
#' should equal the number of suspect sites. If this is present, the function calculates and returns the Gini coefficient 
#' based on the number of crimes; otherwise, this is calculated based on the number of suspect sites.
#'
#' @param hit_scores object in the format defined by geoReportHitscores().
#' @param crimeNumbers optional vector with numbers of crimes per suspect site.
#'
#' @export
#' @examples
#' # John Snow cholera data
#' data(Cholera)
#' d <- geoData(Cholera[,1],Cholera[,2])
#' data(WaterPumps)
#' s <- geoDataSource(WaterPumps[,1], WaterPumps[,2])
#' p <- geoParams(data = d, sigma_mean = 1.0, sigma_squared_shape = 2)
#' m <- geoMCMC(data = d, params = p, lambda=0.05)
#' geoPlotMap(params = p, data = d, source = s, surface = m$geoProfile, breakPercent = seq(0, 50, 5), mapType = "hybrid",
#' crimeCol = "black", crimeCex = 2, sourceCol = "red", sourceCex = 2)
#' hs <- geoReportHitscores(mcmc=m, source=s)
#' # Lorenz plot
#' geoPlotLorenz(hit_scores=hs)
#' 
#' # simulated data
#' sim <-rDPM(50, priorMean_longitude = -0.04217491, priorMean_latitude = 
#' 51.5235505, alpha=1, sigma=1, tau=3)
#' d <- geoData(sim$longitude, sim $latitude)
#' s <- geoDataSource(sim$source_lon, sim$source_lat)
#' p <- geoParams(data = d, sigma_mean = 1.0, sigma_squared_shape = 2)
#' m <- geoMCMC(data = d, params = p)
#' geoPlotMap(params = p, data = d, source = s, surface = m$geoProfile, breakPercent = seq(0, 30, 5), mapType = "terrain", 
#' contourCols=c("blue","white"), crimeCol="black", crimeBorderCol="white", crimeCex=2,
#' sourceCol = "red", sourceCex = 2, opacity = 0.7)
#' hs <- geoReportHitscores(mcmc=m, source=s)
#' # Lorenz plot using number of incidents per source
#' cr <- table(sim$group)
#' geoPlotLorenz(hit_scores=hs,crimeNumbers=cr)

geoPlotLorenz <- function(hit_scores, crimeNumbers=NULL, suspects_col="red", crimes_col="blue") {
   		# define function using trapezoid rule
   		tpzd <- function(x,y)
			{
				idx = 2:length(x)
   				return (as.double( (x[idx] - x[idx-1]) %*% (y[idx] + y[idx-1])) / 2)
   			}

	if(is.null(crimeNumbers))
    {
        # sort hit scores
        capture.output(ordered_hs <- hit_scores[order(hit_scores[,3]),][,3])
        
        # cumulative crime sites
        cum_suspect_sites <- seq(1,length(ordered_hs))/length(ordered_hs)
        
        # add initial 0 and terminal 1 to each so plots will begin at (0,0) and end at (1,1) (these are dropped later)
        ordered_hs <- c(0,ordered_hs,1)
        cum_suspect_sites <- c(0,cum_suspect_sites,1)
        
        # calculate gini coefficient
 		auc <- (0.5-tpzd(cum_suspect_sites, ordered_hs))/0.5
 		G <- round(auc,3)        

        # plot
        plot(ordered_hs, cum_suspect_sites,type="n",xlim=c(0,1),ylim=c(0,1),xlab="hit score",ylab="proportion of sources")
        points(ordered_hs, cum_suspect_sites,type="l",col=suspects_col)
        abline(h=seq(0,1,0.2),col="lightgray",lwd=0.4)
        abline(v=seq(0,1,0.2),col="lightgray",lwd=0.4)
        abline(0,1,col="gray")
        text(0.8,0.1,paste0("G (sources) = ",G),cex=0.8)
        
        # output
        return(gini_coefficient=G)
        
    } else
    
    {
        hit_scores <- cbind(hit_scores, crimeNumbers)
        
        # sort hit scores
        capture.output(ordered_hs <- hit_scores[order(hit_scores[,3]),][,3])
        
        # cumulative crime sites
        cum_suspect_sites <- seq(1,length(ordered_hs))/length(ordered_hs)
        
        # cumulative crime numbers
        cum_crimes <- cumsum(hit_scores[,4])/sum(hit_scores[,4])
        
        # add initial 0 and terminal 1 to each so plots will begin at (0,0) and end at (1,1) (these are dropped later)
        ordered_hs <- c(0,ordered_hs,1)
        cum_crimes <- c(0,cum_crimes,1)
        cum_suspect_sites <- c(0,cum_suspect_sites,1)
        
        # calculate gini coefficient
		auc <- (0.5-tpzd(cum_crimes, ordered_hs))/0.5
 		G <- round(auc,3)         

# plot
        plot(ordered_hs,cum_suspect_sites,type="n",xlim=c(0,1),ylim=c(0,1),xlab="hit score",ylab="proportion of sources and incidents")
        abline(h=seq(0,1,0.2),col="lightgray",lwd=0.4)
        abline(v=seq(0,1,0.2),col="lightgray",lwd=0.4)
        points(ordered_hs, cum_suspect_sites,type="l",col= suspects_col)
        points(ordered_hs, cum_crimes,type="l",col= crimes_col)

        abline(0,1,col="darkgray")
        text(0.85,0.3,paste0("G (crimes) = ",G),cex=0.8)
		legend(0.7,0.2,c("sources","incidents"),col=c(suspects_col,crimes_col),lwd=1,cex=0.8)
        
        # output
        return(gini_coefficient=G)

    }
    
}


#------------------------------------------------
#' Calculate and plot probability of coallocation
#'
#' Calculates the probability that two crimes are from the same source.
#' Also allows an optional plot showing the probabilities of allocation different
#' sources for each of the two selected crimes, using myMCMC$allocation.
#' The function returns the position of the two chosen crimes in the original
#'list, their lon/lat and the probability that they come from the same source.
#'
#' @param crime1 numerical index of first crime.
#' @param crime2 numerical index of second crime.
#' @param coallocation_matrix matrix of coallocations between all observations, as produced the "allocation" output of the function geoMCMC().
#' @param offset vertical offset of second line to ensure readability.
#' @param plot.graph whether to plot the graph (if FALSE simply prints probability to console).
#'
#' @export
#' @examples
#' # John Snow cholera data
#' data(Cholera)
#' d <- geoData(Cholera[,1],Cholera[,2])
#' data(WaterPumps)
#' s <- geoDataSource(WaterPumps[,1], WaterPumps[,2])
#' p <- geoParams(data = d, sigma_mean = 1.0, sigma_squared_shape = 2)
#' m <- geoMCMC(data = d, params = p, lambda=0.05)
#' prob_coallocation(crime1=1, crime2=3, coallocation_matrix=m$allocation)
#' 
#' # simulated data
#' sim <-rDPM(50, priorMean_longitude = -0.04217491, priorMean_latitude = 
#' 51.5235505, alpha=1, sigma=1, tau=3)
#' d <- geoData(sim$longitude, sim $latitude)
#' s <- geoDataSource(sim$source_lon, sim$source_lat)
#' p <- geoParams(data = d, sigma_mean = 1.0, sigma_squared_shape = 2)
#' m <- geoMCMC(data = d, params = p)
#' prob_coallocation(crime1=1, crime2=25, coallocation_matrix=m$allocation)

prob_coallocation <- function(crime1, crime2, coallocation_matrix, offset=0.005, plot.graph=TRUE)
	{
        # convert input to matrix and get dimension
        cmat <- as.matrix(coallocation_matrix)
        d <- ncol(cmat)

        # calculate probability that two crimes come from the same source
        # set to 1 if we are considering the probability that a source coallocates with itself!
        ifelse(crime1==crime2,
        prob_coall<-1,
        prob_coall <- sum(cmat[crime1,]*cmat[crime2,]))

		# plot graph if required
		if(plot.graph)
			{
				# plot lines
                mainTitle <- paste0("comparing crimes ",crime1," and ",crime2)
				plot(cmat[crime1,], type="l", xlim=c(0,d), ylim=c(0,1), col="red", xlab="source", ylab="probability", main=mainTitle)
				points(cmat[crime2,]+offset, type="l", col="darkgray")
                
                # add legend text etc.
				text(d/1.5, 0.9, paste("probability of coallocation =", round(prob_coall,3)), cex=0.9)
				legend(d/1.3,0.2, legend=c(crime1,crime2), lwd=1, col=c("red","darkgray"), cex=0.7)
			}
		
        # return coallocation probability of these two crimes
		return(list("crime.1"=crime1,"crime.2"=crime2,"p.coallocation"=prob_coall))
	}

#------------------------------------------------
#' Calculate and plot hit scores based on a ring search
#'
#' Calculates hit scores for a ring-search strategy (ie searching in an expanding radius out from the crimes). Also plots the crimes and sources with merged polygons showing these (merged and clipped) rings
#'
#' @param params parameters list in the format defined by geoParams().
#' @param data data object in the format defined by geoData().
#' @param source potential sources object in the format defined by geoDataSource().
#' @param buffer_radii vector of search radiuses to draw around incidents, in metres. Default is 1000m, 2000m and 5000.
#'
#' @export
#' @examples
#' # John Snow cholera data
#' data(Cholera)
#' d <- geoData(Cholera[,1],Cholera[,2])
#' data(WaterPumps)
#' s <- geoDataSource(WaterPumps[,1], WaterPumps[,2])
#' p <- geoParams(data = d, sigma_mean = 1.0, sigma_squared_shape = 2)
#' m <- geoMCMC(data = d, params = p, lambda=0.05)
#' geoPlotMap(params = p, data = d, source = s, surface = m$geoProfile, breakPercent = seq(0, 50, 5), mapType = "hybrid",
#' crimeCol = "black", crimeCex = 2, sourceCol = "red", sourceCex = 2)
#' geoReportHitscores(mcmc=m, source=s)
#' ringHS(params = p, data = d, source = s, buffer_radii=c(1000,2000,5000))
#'
#' # simulated data
#' sim <-rDPM(50, priorMean_longitude = -0.04217491, priorMean_latitude = 
#' 51.5235505, alpha=1, sigma=1, tau=3)
#' d <- geoData(sim$longitude, sim $latitude)
#' s <- geoDataSource(sim$source_lon, sim$source_lat)
#' p <- geoParams(data = d, sigma_mean = 1.0, sigma_squared_shape = 2)
#' m <- geoMCMC(data = d, params = p)
#' geoPlotMap(params = p, data = d, source = s, surface = m$geoProfile, breakPercent = seq(0, 30, 5), mapType = "terrain", 
#' contourCols=c("blue","white"), crimeCol="black", crimeBorderCol="white", crimeCex=2,
#' sourceCol = "red", sourceCex = 2, opacity = 0.7)
#' hs <- geoReportHitscores(mcmc=m, source=s)
#' ringHS(params = p, data = d, source = s, buffer_radii=c(1000,2000,5000))

ringHS <- function(params, data, source, buffer_radii=c(1000,2000,5000))
	{
		
		long2UTM <- function(long) {
        	(floor((long + 180)/6)%%60) + 1
		}
		
		my_UTM <- long2UTM(mean(data$longitude))
		my_crs_long_lat <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
		my_crs_utm <- paste("+proj=utm +zone=", my_UTM, " ellps=WGS84", sep = "")
		
		crimes <- cbind(data$longitude, data$latitude)
		crimes_lonlat <- SpatialPointsDataFrame(coords = as.matrix(crimes), data = as.data.frame(crimes), proj4string = CRS(my_crs_long_lat))
		
		sources <- cbind(source$source_longitude, source$source_latitude)
		sources_lonlat <- SpatialPointsDataFrame(coords = as.matrix(sources), data = as.data.frame(sources), proj4string = CRS(my_crs_long_lat))
		
		lonMin <- params$output$longitude_minMax[1]
		lonMax <- params$output$longitude_minMax[2]
		latMin <- params$output$latitude_minMax[1]
		latMax <- params$output$latitude_minMax[2]
		
		bounds = matrix(c(lonMin, latMin, lonMin, latMax, lonMax, latMax, lonMax, latMin, lonMin, latMin), ncol = 2, byrow = TRUE)
		bounds_lonlat <- SpatialPointsDataFrame(coords = as.matrix(bounds), data = as.data.frame(bounds), proj4string = CRS(my_crs_long_lat))
		bounds_UTM <- spTransform(bounds_lonlat, CRS(my_crs_utm))
		p1 = Polygon(bounds_UTM)
		bounds_polygon_utm = SpatialPolygons(list(Polygons(list(p1), ID = "a")), proj4string = CRS(my_crs_utm))
		map_area_m_sq <- gArea(bounds_polygon_utm)
		
		crimes_UTM <- spTransform(crimes_lonlat, CRS(my_crs_utm))
		sources_UTM <- spTransform(sources_lonlat, CRS(my_crs_utm))
		
		no_sources <- dim(sources_lonlat)[1]
		min_dists <- rep(NA, no_sources)
		for (i in 1:no_sources)
			{
				min_dists[i] <- min(spDistsN1(crimes_lonlat, sources_lonlat[i, 1], longlat = TRUE)) * 1000
			}
		stored_results <- matrix(rep(NA, no_sources * 4), ncol = 4)
		stored_buffers <- list()
		colnames(stored_results) <- c("ring_lon", "ring_lat", "merged_buffer_area_m2", "ring_hs")
    
		for (source_number_to_check in 1:no_sources)
			{
				ifelse(
						min_dists[source_number_to_check]==0,
       
						{stored_results[source_number_to_check, ] <- c(sources[source_number_to_check, ], 0, 0)},
       
       					{	b <- gBuffer(crimes_UTM, byid = FALSE, width = min_dists[source_number_to_check])
        					b2 <- gUnaryUnion(b)
        					clip <- gIntersection(b2, bounds_polygon_utm, byid = TRUE, drop_lower_td = TRUE)
							merged_area <- gArea(clip)
							ring_hs_for_this_source <- merged_area/map_area_m_sq
							stored_results[source_number_to_check, ] <- c(sources[source_number_to_check, ], merged_area, ring_hs_for_this_source)
        				}
        				)
    		}
    
		stored_contours <- list()
			for (contour_number in 1:length(buffer_radii))
			{
				b <- gBuffer(crimes_UTM, byid = FALSE, width = buffer_radii[contour_number])
				b2 <- gUnaryUnion(b)
				clip <- gIntersection(b2, bounds_polygon_utm, byid = TRUE, drop_lower_td = TRUE)
				clip_lonlat <- spTransform(clip, CRS(my_crs_long_lat))
				stored_contours[contour_number] <- clip_lonlat
    		}
    
		my_bounds  <- make_bbox(lon = bounds[,1], lat = bounds[,2], f = 1)
		df_data <- data.frame(long=data$longitude, lat=data$latitude)
		source_df_data <- data.frame(long=source$source_longitude, lat=source$source_latitude)
		
		MyMap <- get_map(location=my_bounds)
		MyMap <-ggmap(MyMap)

		gp.colors <- colorRampPalette(c("red","white"))
    	ringCols <- gp.colors(length(stored_contours))
    	ringColsTransp = AddAlpha(ringCols, 0.05)
    
    	for (cc in length(stored_contours):1)
    		{
				MyMap <- MyMap + geom_polygon(aes(x=long, y=lat, group=group), fill=ringColsTransp[cc], color='black', lwd=0.5,data= stored_contours[[cc]],alpha=0.2)
			}
		MyMap  <- MyMap + geom_point(aes(x=long, y=lat),data=df_data) 
		MyMap  <- MyMap + geom_point(aes(x=long, y=lat),data=source_df_data,pch=15,col="blue")

		box_to_plot <- spTransform(bounds_polygon_utm,CRS(my_crs_long_lat))

		MyMap <- MyMap + geom_polygon(aes(x=long, y=lat, group=group),fill=NA, col="black",lwd=0.5,data= box_to_plot,alpha=0.2)

		plot(MyMap)    
      
    ring_table <- data.frame(stored_results)
    ring_hs <- ring_table[, 4]
    ring_areas <- ring_table[, 3]
    ring_output <- list(ring_table = ring_table, ring_areas = ring_areas, ring_hs = ring_hs, my_UTM = my_UTM, map_area_m_sq = map_area_m_sq)
    return(ring_output)
}

#------------------------------------------------
#' Calculate and plot hit scores based on a ring search
#'
#' Second attempt at ring search without using other packages. TODO - complete help for this function!
#'
#' @param params Parameters list in the format defined by geoParams().
#' @param data Data object in the format defined by geoData().
#' @param source Potential sources object in the format defined by geoDataSource().
#' @param mcmc mcmc object of the form produced by geoMCMC(). 
#' @param buffer_radii Optional vector giving diameter of rings (in km) to show around each crime, suitable merged and clipped. If NULL, full contour map of ring hitscores is plotted instead. 
#'
#' @export
#' @examples
#' # John Snow cholera data
#' data(Cholera)
#' d <- geoData(Cholera[,1],Cholera[,2])
#' data(WaterPumps)
#' s <- geoDataSource(WaterPumps[,1], WaterPumps[,2])
#' p <- geoParams(data = d, sigma_mean = 1.0, sigma_squared_shape = 2)
#' m <- geoMCMC(data = d, params = p, lambda=0.05)
#' geoPlotMap(params = p, data = d, source = s, surface = m$geoProfile, breakPercent = seq(0, 50, 5), mapType = "hybrid",
#' crimeCol = "black", crimeCex = 2, sourceCol = "red", sourceCex = 2)
#' geoReportHitscores(mcmc=m, source=s)
#' # plot full contour plot of ring search
#' geoRingHitscores(params = p, data = d, source = s, mcmc = m)
#' # plot rings of 100m, 200m and 500m around points
#' geoRingHitscores(params = p, data = d, source = s, mcmc = m, buffer_radii = c(0.1, 0.2, 0.5))
#'

geoRingHitscores <- function(params, data, source, mcmc, buffer_radii = NULL) {
    
    # Calculates the percentage of the grid that must be searched before reaching each source under a ring search strategy. This search strategy assumes that we start from a given crime and search outwards in a circle of increasing radius until we reach a source. As there are multiple crimes the strategy assumes a separate individual searching from each crime simultaneously at an equal rate.
    # The basic logic of the approach here is that calculating the final search radius is easy - it is simply the minimum radius from any crime to this source. The difficulty is calculating the amount of grid that will have been explored by the time we reach this radius, as circles will often overlap and the intersection should not be double-counted (we assume searching moves on if the area has already been explored by someone else). This is done by brute force - a grid is created and cells are filled in solid if they have been explored. The total percentage of filled cells gives the hitscore percentage. The distance matrices used in this brute force step are needed repeatedly, and so they are computed once at the begninning to save time.
    
    # get number of crimes and sources
    n <- length(data$latitude)
    ns <- length(source$source_latitude)
    
    # create matrices giving lat/lon at all points in search grid
    lonVec <- mcmc$midpoints_longitude
    latVec <- mcmc$midpoints_latitude
    lonMat <- matrix(rep(lonVec,each=length(latVec)), length(latVec))
    latMat <- matrix(rep(latVec,length(lonVec)), length(latVec))
    
    # calculate great circle distance from every data point to every point in search grid. This list of distance matrices will be used multiple times so best to pre-compute here.
    dlist <- list()
    for (i in 1:n) {
        dlist[[i]] <- latlon_to_bearing(data$latitude[i], data$longitude[i], latMat, lonMat)$gc_dist
    }
    
    # calculate hitscore for each source in turn
    hitScore <- rep(NA,ns)
    for (i in 1:ns) {
        
        # get minimum distance (in km) from all crimes to this source. This is how far the ring search must go out before it reaches the source
        minDist <- min(latlon_to_bearing(source$source_latitude[i], source$source_longitude[i], data$latitude, data$longitude)$gc_dist)
        
        # loop through all crimes, filling in a search matrix as we go if cells are within this minimum distance
        searchMat <- matrix(0,length(latVec),length(lonVec))
        for (j in 1:n) {
            searchMat[dlist[[j]]<=minDist] <- 1
        }
        
        # calculate hitscore percentage from proportion of filled cells
        hitScore[i] <- mean(searchMat)
        
    }
    # convert dlist to an array
    darray <- array(as.vector(unlist(dlist)), dim = c(length(lonVec),length(latVec),n))
    nearest_crime_dist <- t(apply(darray,c(1,2),min))
    
    # if buffer_radii is supplied, plot map showing this. If not, plot 'geoprofile' showing ring search
ifelse(is.null(buffer_radii)==TRUE,
   {print(geoPlotMap(data = data, source = source, params = p, breakPercent = seq(0, 100, 20), mapType = "roadmap", contourCols =c("red", "white"), crimeCol = "black", crimeCex = 2, sourceCol = "red", sourceCex = 2, surface = rank(nearest_crime_dist)))},
    {
ringsMat <- matrix(rep(length(buffer_radii)+1,(params$output$longitude_cells*params$output$latitude_cells)),nrow=params$output$latitude_cells)
buffer_radii <- buffer_radii[order(buffer_radii)]
for(b in length(buffer_radii):1)
{for(i in 1:params$output$longitude_cells)
{
	for(j in 1:params$output$latitude_cells)
	{
		if(nearest_crime_dist[i,j] <= buffer_radii[b]) ringsMat[i,j] <- b
	
	}
}
}
# contour(ringsMat)
print(geoPlotMap(data = d, source = s, params = p, breakPercent = seq(0, 100, 5), mapType = "roadmap", contourCols =c("red", "white"), crimeCol = "black", crimeCex = 2, sourceCol = "red", sourceCex = 2, surface = rank(ringsMat),gpLegend=FALSE))})
    
    ringHS <- hitScore
    ringResults <- data.frame(cbind(source$source_longitude,source$source_latitude, ringHS))
    colnames(ringResults) <- c("lon","lat","ringHS")    
    return(ringResults)
}

#------------------------------------------------
#' Perspective plot of geoprofile or raw probabilities
#'
#' Plots persp plot of geoprofile or posterior surface (coloured according to height), reducing matrix dimensions if necessary to avoid grid lines being too close together. NB Only works with square matrix
#'
#' @param surface surface to plot; either the geoprofile or posteriorSurface output by geoMCMC(). 
#' @param aggregate_size the number of cells to aggregate to smooth the surface.
#' @param surface_type type of surface; should be either "gp" for geoprofile or "prob" for posteriorSurface.
#' @param perspCol colour palette. Defaults to red/orange/yellow/white.
#' @param phiGP value of phi to pass to persp().
#' @param thetaGP value of theta to pass to persp().
#'
#' @export
#' @examples
#' # John Snow cholera data
#' data(Cholera)
#' d <- geoData(Cholera[,1],Cholera[,2])
#' data(WaterPumps)
#' s <- geoDataSource(WaterPumps[,1], WaterPumps[,2])
#' p <- geoParams(data = d, sigma_mean = 1.0, sigma_squared_shape = 2)
#' m <- geoMCMC(data = d, params = p, lambda=0.05)
#' # raw probabilities
#' perspGP(m$posteriorSurface, surface_type = "prob")
#' # geoprofile
#' perspGP(m$geoProfile, aggregate_size = 3, surface_type = "gp")
#' 
#' # simulated data
#' sim <-rDPM(50, priorMean_longitude = -0.04217491, priorMean_latitude = 
#' 51.5235505, alpha=1, sigma=1, tau=3)
#' d <- geoData(sim$longitude, sim $latitude)
#' s <- geoDataSource(sim$source_lon, sim$source_lat)
#' p <- geoParams(data = d, sigma_mean = 1.0, sigma_squared_shape = 2)
#' m <- geoMCMC(data = d, params = p)
#' # raw probabilities
#' perspGP(m$posteriorSurface, surface_type = "prob")
#' # geoprofile
#' perspGP(surface = m$geoProfile, aggregate_size = 3, surface_type = "gp")

perspGP <- function(surface, aggregate_size = 3, surface_type = "gp", perspCol = c("red", "orange", "yellow", "white"), phiGP = 30, thetaGP = -30) {
    
			matrix_manipulation <- function(my_matrix,my_operation)
				{
					if(my_operation=="reflect_y") {return(my_matrix[,ncol(my_matrix):1])}
					if(my_operation=="reflect_x") {return(my_matrix[nrow(my_matrix):1,])}
					if(my_operation=="rotate_180") {return(my_matrix[nrow(my_matrix):1,ncol(my_matrix):1])}
					if(my_operation=="reflect_diag") {return(t(my_matrix))}
					if(my_operation=="rotate_90") {return(t(my_matrix)[,nrow(my_matrix):1])}
					if(my_operation=="rotate_-90") {return(t(my_matrix)[ncol(my_matrix):1,])}
				}

			if(surface_type=="gp") scale <- -1
			if(surface_type=="prob") scale <- 1
			
			
			
			to_plot <-  matrix_manipulation(scale*surface,"reflect_diag")
			
			# reduce matrix or not
			matrix_size <- unique(dim(to_plot))[1]
			breaks <- seq(1,(matrix_size-(aggregate_size-1)), aggregate_size)
			output <- matrix(rep(NA,length(breaks)^2),ncol=length(breaks))
			for(i in 1: length(breaks))
				{
					for(j in 1: length(breaks))
						{
							output[i,j] <- mean(as.vector(to_plot[breaks[i]:(breaks[i]+(aggregate_size-1)),breaks[j]:(breaks[j]+(aggregate_size-1))]))
						}
				}
                
			# select colours
			gp.colors <- colorRampPalette(perspCol)
            
			# Generate the desired number of colors from this palette
			nbcol <- 100
			color <- gp.colors(nbcol)
			ncz <- dim(output)[2]
			nrz <- dim(output)[1]
            
			# Compute the z-value at the facet centres
			zfacet <- output[-1, -1] + output[-1, -ncz] + output[-nrz, -1] + output[-nrz, -ncz]
			facetcol <- cut(zfacet, nbcol)
		
			persp(output, col=color[facetcol], border="black", phi=phiGP, theta=thetaGP, lwd=0.2, box=FALSE)
            
	}


#------------------------------------------------
#' Interactive perspective plot of geoprofile or raw probabilities
#'
#' Produces interactive perspective plot of geoprofile or posterior surface (coloured according to height).
#'
#' @param surface surface to plot; either the geoprofile or posteriorSurface output by geoMCMC().
#' @param surface_type type of surface; should be either "gp" for geoprofile or "prob" for posteriorSurface.
#' @param perspCol colour palette. Defaults to red/orange/yellow/white.
#' @param scale vertical scale of surface.
#'
#' @export
#' @examples
#' # John Snow cholera data
#' data(Cholera)
#' d <- geoData(Cholera[,1],Cholera[,2])
#' data(WaterPumps)
#' s <- geoDataSource(WaterPumps[,1], WaterPumps[,2])
#' p <- geoParams(data = d, sigma_mean = 1.0, sigma_squared_shape = 2)
#' m <- geoMCMC(data = d, params = p, lambda=0.05)
#' # raw probabilities
#' perspGP2(m$posteriorSurface, surface_type = "prob")
#' # geoprofile
#' perspGP2(m$geoProfile, surface_type = "gp")
#' 
#' # simulated data
#' sim <-rDPM(50, priorMean_longitude = -0.04217491, priorMean_latitude = 
#' 51.5235505, alpha=1, sigma=1, tau=3)
#' d <- geoData(sim$longitude, sim $latitude)
#' s <- geoDataSource(sim$source_lon, sim$source_lat)
#' p <- geoParams(data = d, sigma_mean = 1.0, sigma_squared_shape = 2)
#' m <- geoMCMC(data = d, params = p)
#' # raw probabilities
#' perspGP2(m$posteriorSurface, surface_type = "prob")
#' # geoprofile
#' perspGP2(surface = m$geoProfile, surface_type = "gp")

perspGP2 <- function(surface, surface_type="gp", perspCol=c("red", "orange", "yellow", "white"), scale=1) {
    
    # check input formats
    stopifnot(surface_type%in%c("gp","prob"))
    
    # produce x, y and z values
    xvec <- seq(1,0,l=nrow(surface))
    yvec <- seq(0,1,l=ncol(surface))
    zmat <- (surface-min(surface))/(max(surface)-min(surface))
    
    # invert if using gp surface
    if (surface_type=="gp") {
        zmat <- 1-zmat
    }
    
    # produce colors
    ncols <- 1e2
    colRamp <- colorRampPalette(perspCol)
    colMat <- colRamp(ncols)[cut(zmat,ncols)]
    
    # plot surface
    surface3d(x=xvec, y=yvec, z=scale*zmat, color=colMat)
    
}

#------------------------------------------------
#' Incorporate shapefile or raster information into a geoprofile
#' 
#' This function allows information from a shapefile or raster to be incorporated within the geoprofile. For example, we might wish to exclude areas not on land, or weight the probabilities within a specific postcode differently. The spatial object used should be a SpatialPolygonsDataFrame as produced by the package sp or a raster. 
#' 
#' @param probSurface the original geoprofile, usually the object $posteriorSurface produced by geoMCMC().
#' @param params an object produced by geoParams().
#' @param shapefile the spatial information to include. Must be either SpatialPolygonsDataFrame or RasterLayer.
#' @param masterProj the projection to use, eg "+proj=longlat +datum=WGS84".
#' @param scaleValue different functions depending on value of "operation". For "inside' or "outside", the value by which probabilities should be multiplied inside or outside the shapefile; set to zero to exclude completely. For "near" and "far", the importance of proximity to, or distance from, the object described in the RasterLayer or SpatialPointsDataFrame. Not used for "continuous".
#' @param operation thow to combine the surface and the new spatial information. Must be one of "inside", "outside", "near", "far" or "continuous". The first two multiply areas inside or outside the area described in the shapefile (or raster) by scaleValue. "near" or "far" weight the geoprofile by its closeness to (or distance from) the area described in the shapefile (or raster). Finally, "continuous" uses a set of numerical values (eg altitude) to weight the geoprofile.
#' @param maths one of "add", "subtract", multiply" or "divide. The mathematical operation used to combine the new spatial data with the geoprofile when operation = "continuous".
#' 
#' @export
#' @examples
#' # to come

GPshapefile <- function (probSurface, params, shapefile, masterProj = "+proj=longlat +datum=WGS84", 
    scaleValue = 1, operation = "inside", maths = "multiply") 
{
    stopifnot(class(shapefile) %in% c("SpatialPolygonsDataFrame", 
        "RasterLayer"))
    stopifnot(operation %in% c("inside", "outside", "near", "far", 
        "continuous"))
    stopifnot(maths %in% c("multiply", "divide", "add", "subtract", 
        "continuous"))
 if (class(shapefile) == "RasterLayer") {
       rf <- shapefile
    }
    if (class(shapefile) == "SpatialPolygonsDataFrame") {
        r <- raster(ncol = params$output$longitude_cells, nrow = params$output$latitude_cells)
    extent(r) <- extent(shapefile)
    rf <- rasterize(shapefile, r)
     }
probSurface <- probSurface[params$output$longitude_cells:1, 
        ]
    master_extent <- rbind(params$output$longitude_minMax, params$output$latitude_minMax)
    masterproj <- masterProj
    r <- raster(probSurface, crs = masterproj, xmn = params$output$longitude_minMax[1], 
        xmx = params$output$longitude_minMax[2], ymn = params$output$latitude_minMax[1], 
        ymx = params$output$latitude_minMax[2])
    raster_probSurface <- r
 new_spatial_data_to_include <- projectRaster(rf, raster_probSurface, 
        crs = master_proj)
    new_data_as_scaled_matrix <- matrix(new_spatial_data_to_include@data@values, 
        ncol = params$output$longitude_cells)[, params$output$latitude_cells:1]
    GP_as_scaled_matrix <- matrix(raster_probSurface@data@values, 
        ncol = params$output$longitude_cells)[, params$output$latitude_cells:1]
    combined_mat <- matrix(rep(NA, (params$output$longitude_cells * 
        params$output$latitude_cells), nrows = params$output$longitude_cells), 
        ncol = params$output$longitude_cells)
 if (operation == "inside") {
        for (i in 1:params$output$longitude_cells) {
            for (j in 1:params$output$latitude_cells) {
                new_value <- ifelse(is.na(new_data_as_scaled_matrix[i, 
                  j]) == FALSE, GP_as_scaled_matrix[i, j], (scaleValue * 
                  GP_as_scaled_matrix[i, j]))
                combined_mat[i, j] <- new_value
            }
        }
        scaleMatrix <- new_data_as_scaled_matrix
    }
    if (operation == "outside") {
        for (i in 1:params$output$longitude_cells) {
            for (j in 1:params$output$latitude_cells) {
                new_value <- ifelse(is.na(new_data_as_scaled_matrix[i, 
                  j]) == TRUE, GP_as_scaled_matrix[i, j], (scaleValue * 
                  GP_as_scaled_matrix[i, j]))
                combined_mat[i, j] <- new_value
            }
        }
        scaleMatrix <- new_data_as_scaled_matrix
    }
    if (operation == "continuous") {
        if (maths == "add") {
            combined_mat <- GP_as_scaled_matrix + new_data_as_scaled_matrix
        }
        if (maths == "subtract") {
            combined_mat <- GP_as_scaled_matrix - new_data_as_scaled_matrix
        }
        if (maths == "multiply") {
            combined_mat <- GP_as_scaled_matrix * new_data_as_scaled_matrix
        }
        if (maths == "divide") {
            combined_mat <- GP_as_scaled_matrix/new_data_as_scaled_matrix
        }
        scaleMatrix <- new_data_as_scaled_matrix
    }
    if (operation == "near") {
        distance_mat <- distance(new_spatial_data_to_include)
        distance_mat <- matrix(distance_mat@data@values, ncol = params$output$longitude_cells)[, 
            params$output$latitude_cells:1]
        scaleMatrix <- 1/(distance_mat^scaleValue)
        scaleMatrix[scaleMatrix == "Inf"] <- 1
        combined_mat <- scaleMatrix * GP_as_scaled_matrix
    }
    if (operation == "far") {
        distance_mat <- distance(new_spatial_data_to_include)
        distance_mat <- matrix(distance_mat@data@values, ncol = params$output$longitude_cells)[, 
            params$output$latitude_cells:1]
        scaleMatrix <- distance_mat^scaleValue
        scaleMatrix[scaleMatrix == "Inf"] <- 1
        combined_mat <- scaleMatrix * GP_as_scaled_matrix
    }
    adjusted_surface <- combined_mat
    rank_adjusted_surface <- rank(-adjusted_surface)
    adjSurface <- list(rank = matrix(rank_adjusted_surface, ncol = params$output$longitude_cells, 
        byrow = TRUE), prob = matrix(adjusted_surface, ncol = params$output$longitude_cells, 
        byrow = TRUE)/sum(adjusted_surface[is.na(adjusted_surface)==FALSE]), scaleMatrix = scaleMatrix)
    return(adjSurface)
   }
#------------------------------------------------
#' Extract latitude and longitude of points identified as sources by geoMCMC()
#' 
#' This function takes the output of geoMCMC() and, for each 'crime', extracts the group to which it is assigned with the highest probability. For each group, the model returns a list of these groups, and the mean lat/long of all crimes assigned to that group, returning these in a format matching the output of geoDataSource() for easy plotting with geoPLotMap().
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
#' ms <- modelSources(mcmc = m, data = d)
#' # plot data showing the sources identified by the model (note: NOT the actual suspect sites)
#' geoPlotMap(data = d, source = ms$model_sources, params = p, breakPercent = seq(0, 10, 1), 
#' mapType = "roadmap", contourCols =c("red", "orange","yellow","white"), crimeCol = "black",
#' crimeCex = 2, sourceCol = "red", sourceCex = 2, surface = m$geoProfile, gpLegend=TRUE,
#' opacity = 0.4)

modelSources <- function (mcmc, data) 
{
    groups <- apply(mcmc$allocation, 1, which.max)
    group_IDs <- unique(groups)
    max_groups <- max(group_IDs)
    sources_found <- data.frame(matrix(rep(NA, 2 * max_groups), ncol = 2))
    n_groups <- length(group_IDs)
    for (i in 1:n_groups)
    {
        sources_found[group_IDs[i], ] <- c(mean(data$longitude[which(groups == group_IDs[i])]), mean(data$latitude[which(groups == group_IDs[i])]))
    }
    sources_found <- sources_found[complete.cases(sources_found),]
    model_sources <- geoDataSource(sources_found[, 1],sources_found[,2])
    return(list(model_sources = model_sources, groups = groups))
}
#------------------------------------------------
#' Interactive map with zoom function
#' 
#' Like geoPlotMap(), this function takes the output of geoMCMC() and plots the resulting geoprofile, but this time in an active window allowing the user to click on two points that define an area on which to zoom. For simplicity, the original plot and the zoom plot returned by the function lack most of the options of geoPlotMap() (for example, it doesn't plot sources, or allow custom colours, contours etc). However, the function returns new params and surface objects ($paramsZoom and $surfaceZoom respectively) which can be used with geoPLotMap(). NOTE: The function gglocator() from ggmap is relatively slow to process the first click, so users should wait until the cross hairs reappear before clicking a second time.
#' 
#' @param my_data Crime site data, in the format produced by geoData().
#' @param my_params A params object, in the format produced by geoParams().
#' @param my_surface Surface to plot.
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
#' z <- geoPlotZoom(my_data = d, my_params = p, my_surface = m$geoProfile)
#' # replot zoom, customising the output with geoPlotMap()
#' geoPlotMap(data = d, source = s, params = z$paramsZoom, breakPercent = seq(0, 10, 1), 
#' mapType = "roadmap", contourCols =c("red", "orange","yellow","white"), 
#' crimeCol = "black", crimeCex = 2, sourceCol = "red", sourceCex = 2, 
#' surface = z$surfaceZoom, gpLegend=TRUE)

geoPlotZoom <- function(my_data, my_params, my_surface)
{
# plot original map
print(geoPlotMap(data = my_data, params = my_params, surface = my_surface))

# extract lat long limits for zoom
print("Please click on the active map window to select two points defining the area you wish to magnify:",quote=FALSE)
pos <- gglocator(2)

# create new params object from original params
p2<-p
p2$output$longitude_minMax <- c(min(pos[,1]),max(pos[,1]))
p2$output$latitude_minMax <- c(min(pos[,2]),max(pos[,2]))

# calculate part of original geoprofile to extract
lon_vec_orig <- seq(p$output$longitude_minMax[1],p$output$longitude_minMax[2],length=p$output$longitude_cells)
lat_vec_orig <- seq(p$output$latitude_minMax[1],p$output$latitude_minMax[2],length=p$output$latitude_cells)
mat_min_lon <- which(lon_vec_orig>p2$output$longitude_minMax[1])[1]
mat_max_lon <- which(lon_vec_orig>p2$output$longitude_minMax[2])[1]
mat_min_lat <- which(lat_vec_orig>p2$output$latitude_minMax[1])[1]
mat_max_lat <- which(lat_vec_orig>p2$output$latitude_minMax[2])[1]
sub_mat <- my_surface[mat_min_lat:mat_max_lat,mat_min_lon: mat_max_lon]

# function to resize this sub-matrix to the original resolution
expandMatrix <- function(mat,output_long,output_lat)
{
	# define function expanding vector
	expandVector <- function(input_vec,output_length)
		{
			my_vec <- input_vec
			desired_length <- output_length
			new_vec <- rep(NA, desired_length)

			vec_ID <- seq(1,length(my_vec),length=desired_length)

			for(i in 1:length(new_vec))
				{
					ifelse(vec_ID[i] %% 1 == 0,
		new_vec[i] <- my_vec[floor(vec_ID[i])],
		new_vec[i] <- mean((1-vec_ID[i] %% 1) * my_vec[floor(vec_ID[i])] + (vec_ID[i] %% 1) * my_vec[ceiling(vec_ID[i])])
	)
				}
return(new_vec)
		}
mat1 <- apply(mat,2, function(x) expandVector(x, output_long))
mat2 <- apply(mat1,1, function(x) expandVector(x, output_lat))
return(t(mat2))
}

# resize zoom area of geoprofile
zoomed <- expandMatrix(sub_mat,p$output$longitude_cells,p$output$latitude_cells)

# plot without sources
print(geoPlotMap(data = d, params = p2, surface = zoomed))

# return params and surface objects for further plotting if required
return(list(paramsZoom=p2,surfaceZoom=zoomed))
}
#------------------------------------------------
#' Unknown pleasures
#' 
#' A frivolous alternative to geoPlotMap(), this function takes the output of geoMCMC() and plots the resulting geoprofile in the style of the cover of Joy Division's 'Unknown pleasures' album.
#' 
#' @param input_matrix The surface to plot, usually the object $geoProfile produced by geoMCMC().
#' @param nlines The number of lines (defaults to the correct number of 80).
#' @param paper_ref A text string, for example a reference to a paper.
#' @param bgcol Background colour
#' @param fgcol Foreground colour
#' @param wt line weight
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
#' unknownPleasures(m$geoProfile, paper_ref = "Rgeoprofile v2.0.0")

unknownPleasures <- function(input_matrix, paper_ref = NULL, nlines = 80, bgcol = "black", fgcol = "white", wt = 2)
{
	orig_prof <- input_matrix
	citation <- paper_ref
	# extract other params
	ncols <- ncol(orig_prof)
# functions
expandMatrix <- function(mat,output_long,output_lat)
{
	# define function expanding vector
	expandVector <- function(input_vec,output_length)
		{
			my_vec <- input_vec
			desired_length <- output_length
			new_vec <- rep(NA, desired_length)

			vec_ID <- seq(1,length(my_vec),length=desired_length)

			for(i in 1:length(new_vec))
				{
					ifelse(vec_ID[i] %% 1 == 0,
		new_vec[i] <- my_vec[floor(vec_ID[i])],
		new_vec[i] <- mean((1-vec_ID[i] %% 1) * my_vec[floor(vec_ID[i])] + (vec_ID[i] %% 1) * my_vec[ceiling(vec_ID[i])])
	)
				}
return(new_vec)
		}
mat1 <- apply(mat,2, function(x) expandVector(x, output_long))
mat2 <- apply(mat1,1, function(x) expandVector(x, output_lat))
return(t(mat2))	
}
# reduce to manageable number of rows and columns!
reduced_mat <- expandMatrix(orig_prof ,nlines, ncols)
# scale so values fall between -0.5 and +0.5
reduced_mat <- 1-reduced_mat/max(reduced_mat)-0.5
# set y coordinates of lines
yvals <- seq(0, 1,length = nlines)
# plot
par(bg = bgcol)
plot(1:ncols,reduced_mat[nlines,]+yvals[nlines],type="l",ylim=c(-1,max(yvals)+0.5),axes=FALSE,xlab="",ylab="",col=fgcol)
for(i in nlines:1)
{
	
	polygon(c(min(reduced_mat[i,]),reduced_mat[i,]+yvals[i],min(reduced_mat[i,])),col=bgcol,border=bgcol)
	points(1:ncol(reduced_mat),reduced_mat[i,]+yvals[i],type="l", col=fgcol,lwd=wt)
}
text(0,-0.6,citation,adj=0,col=fgcol)
}
#------------------------------------------------
#' Plot a map and overlay data and/or geoprofile
#'
#' Plots geoprofile on map, with various customisable options.
#'
#' @param params parameters list in the format defined by geoParams().
#' @param data data object in the format defined by geoData().
#' @param source potential sources object in the format defined by geoDataSource().
#' @param crimeNames text labels for crimes. If NULL, numbers will be used.
#' @param sourceNames text labels for crimes. If NULL, numbers will be used.
#' @param surface a surface to overlay onto the map, typically a geoprofile obtained from the output of geoMCMC().
#' @param zoom zoom level of map. If NULL then choose optimal zoom from params.
#' @param mapSource which online source to use when downloading the map. Options include Google Maps ("google"), OpenStreetMap ("osm"), Stamen Maps ("stamen") and CloudMade maps ("cloudmade").
#' @param mapType the specific type of map to plot. Options available are "terrain", "terrain-background", "satellite", "roadmap" and "hybrid" (google maps), "terrain", "watercolor" and "toner" (stamen maps) or a positive integer for cloudmade maps (see ?get_cloudmademap from the package ggmap for details).
#' @param opacity value between 0 and 1 givin the opacity of surface colours.
#' @param plotContours whether or not to add contours to the surface plot.
#' @param breakPercent vector of values between 0 and 100 describing where in the surface contours appear.
#' @param contourCols list of two or more colours from which to derive the contour colours.
#' @param crimeCex relative size of symbols showing crimes.
#' @param crimeCol colour of crime symbols.
#' @param crimeBorderCol border colour of crime symbols.
#' @param crimeBorderWidth width of border of crime symbols.
#' @param sourceCex relative size of symbols showing suspect sites.
#' @param sourceCol colour of suspect sites symbols.
#' @param gpLegend whether or not to add legend to plot.
#'
#' @export
#' @examples
#' # John Snow cholera data
#' data(Cholera)
#' d <- geoData(Cholera[,1],Cholera[,2])
#' data(WaterPumps)
#' s <- geoDataSource(WaterPumps[,1], WaterPumps[,2])
#' p <- geoParams(data = d, sigma_mean = 1.0, sigma_squared_shape = 2)
#' m <- geoMCMC(data = d, params = p, lambda=0.05)
#' # produce simple map
#' geoPlotMapText(params = p, data = d, source = s, crimeNames = NULL, sourceNames =
#' letters[1:13], surface = m$geoProfile, breakPercent = seq(0, 50, 5), mapType = "hybrid",
#' crimeCol = "black", crimeCex = 2, sourceCol = "red", sourceCex = 2)
#' 

geoPlotMapText <- function (params, data = NULL, source = NULL, crimeNames = NULL, sourceNames = NULL, surface = NULL, 
    zoom = NULL, mapSource = "google", mapType = "hybrid", opacity = 0.6, 
    plotContours = TRUE, breakPercent = seq(0, 100, l = 11), 
    contourCols = c("red", "orange", "yellow", "white"), crimeCex = 1.5, 
    crimeCol = "red", crimeBorderCol = "white", crimeBorderWidth = 0.5, 
    sourceCex = 1.5, sourceCol = "blue", gpLegend = TRUE) 
    {
    geoParamsCheck(params)
    if (!is.null(data)) 
        geoDataCheck(data)
    if (is.null(zoom)) 
        zoom <- getZoom(params$output$longitude_minMax, params$output$latitude_minMax)
    if (mapSource == "stamen") 
        zoom <- min(zoom, 18)
    rawMap <- get_map(location = c(mean(params$output$longitude_minMax), 
        mean(params$output$latitude_minMax)), zoom = zoom, source = mapSource, 
        maptype = mapType)
    myMap <- ggmap(rawMap) + coord_cartesian(xlim = params$output$longitude_minMax, 
        ylim = params$output$latitude_minMax)
    if (!is.null(surface)) {
        geoCols <- colorRampPalette(contourCols)
        nbcol = length(breakPercent) - 1
        longitude_minMax <- params$output$longitude_minMax
        latitude_minMax <- params$output$latitude_minMax
        longitude_cells <- params$output$longitude_cells
        latitude_cells <- params$output$latitude_cells
        longitude_cellSize <- diff(longitude_minMax)/longitude_cells
        latitude_cellSize <- diff(latitude_minMax)/latitude_cells
        longitude_midpoints <- longitude_minMax[1] - longitude_cellSize/2 + 
            (1:longitude_cells) * longitude_cellSize
        latitude_midpoints <- latitude_minMax[1] - latitude_cellSize/2 + 
            (1:latitude_cells) * latitude_cellSize
        df <- expand.grid(x = longitude_midpoints, y = latitude_midpoints)
        df$z <- as.vector(t(surface))
        labs <- paste(round(breakPercent, 1)[-length(breakPercent)], 
            "-", round(breakPercent, 1)[-1], "%", sep = "")
        df$cut <- cut(df$z, breakPercent/100 * length(surface), 
            labels = labs)
        df_noNA <- df[!is.na(df$cut), ]
        myMap <- myMap + geom_tile(aes(x = x, y = y, fill = cut), 
            alpha = opacity, data = df_noNA)
        myMap <- myMap + scale_fill_manual(name = "Hitscore\npercentage", 
            values = rev(geoCols(nbcol)))
        if (gpLegend == FALSE) {
            myMap <- myMap + theme(legend.position = "none")
        }
        if (plotContours) {
            myMap <- myMap + stat_contour(aes(x = x, y = y, z = z), 
                colour = "grey50", breaks = breakPercent/100 * 
                  length(surface), size = 0.3, alpha = opacity, 
                data = df)
        }
    }
    	if (is.null(crimeNames)) {crimeNames = 1:length(data$longitude)}
    if (!is.null(data)) {
        df_data <- data.frame(longitude = data$longitude, latitude = data$latitude,ptno=crimeNames)
        df_data$ptno <- crimeNames
        myMap <- myMap + geom_text(aes(x = longitude, y = latitude, label = ptno), 
            data = df_data, cex = crimeCex,col=crimeCol)
    }
    	if (is.null(sourceNames)) {sourceNames = 1:length(source$source_longitude)}
      if (!is.null(source)) {
        df_source <- data.frame(longitude = source$source_longitude, 
            latitude = source$source_latitude,sourceNames)
        df_source$sptno <- sourceNames
        myMap <- myMap + geom_text(aes(x = longitude, y = latitude, label = sptno), 
            data = df_source, cex = sourceCex, col = sourceCol)  
            
    }
    myMap
}
#------------------------------------------------

