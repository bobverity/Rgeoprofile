
# TODO
# - the function ringHS() uses a lot of packages and seems quite complicated! At the moment all we did was change the unput variable names, but not the names in the function, meaning it won't currently work. If you could have a go at overhauling this function with the minimum dependencies that would be great. For example, if we can get away with adding our own circles rather than using spatialDataFrame objects that would be ideal.
# - expand help text and examples where needed. Remember that you need to run the function document() to actually create the help files once you've updated the text here. (NB, document() is part of devtools).
# - add Gini coefficient calculation to plotLorenz() function. Use simple trapezoidal rule rather than a package. For example, if you have vectors x and y of same length then the area under curve is given by sum(0.5*(y[-1]+y[-length(y)])*(x[-1]-x[-length(x)]))
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
#' rDPM(10)

# The following commands are needed to ensure that the roxygen2 package, which deals with documenting the package, does not conflict with the Rcpp package. Do not alter!
#' @useDynLib RgeoProfile
#' @importFrom Rcpp evalCpp
#' @import fftwtools
#' @import ggplot2
#' @import ggmap
#' @import RColorBrewer
#' @exportPattern "^[[:alpha:]]+"

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
#' Draw from bivariate normal distribution transformed to spatial coordinates
#'
#' Draw from normal distribution converted to spherical coordinate system. Points are first drawn from an ordinary cartesian 2D normal distribution. The distances to points are then assumed to be great circle distances, and are combined with a random bearing from the point {centre_lat, centre_lon} to produce a final set of lat/lon points. Note that this is not a truly spherical normal distribution, as the domain of the distribution is not the sphere - rather it is a transformation from one coordinate system to another that is satisfactory when the curvature of the sphere is not severe.
#'
#' @param n number of draws.
#' @param centre_lat latitude of the centre point of the normal distribution.
#' @param centre_lon longitude of the centre point of the normal distribution.
#' @param sigma standard deviation of normal distribution in km.
#'
#' @export
#' @examples
#' rnorm_sphere(5, 51.5074, -0.1277, 1)

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
#' @param longitude the longitudinal positions of the observed data.
#' @param latitude the latitudinal positions of the observed data.
#'
#' @export
#' @examples
#' geoData()

geoData <- function(longitude=NULL, latitude=NULL) {
    
  # generate dummy data if none read in
  if (is.null(longitude) & is.null(latitude)) {
    longitude <- c(-0.104545976281589, -0.102659272660916, -0.0967390020136406, -0.0996246226730725, -0.100775342233937, -0.101073477576196, -0.100932674617746, -0.0983001766339886, -0.0913571765598557, -0.100211479242536, -0.139508969429415, -0.14403082311245, -0.143607222414313, -0.137174795971723, -0.140884394738737, -0.142723755125487, -0.143380928147727, -0.136989691342132, -0.13837666855334, -0.138297288871952, -0.0773858357074935, -0.0818743917621333, -0.0738310357273188, -0.0744118149244568, -0.0757833597110897, -0.0762193916493531, -0.0810467015747727, -0.110052994420826, -0.106600836874167, -0.105104028808356, -0.101934241194567, -0.0683752111183375, -0.0758607240702608, -0.079153744918552, -0.087964365345432)
    latitude <- c(51.4996147329979, 51.4925230579844, 51.4947129689414, 51.4922683109254, 51.5007532206834, 51.4960640374896, 51.4996917836745, 51.4976936749008, 51.4977904998888, 51.4894186202378, 51.5002583182117, 51.5033510595094, 51.4984697991335, 51.5063306206839, 51.4961516950408, 51.4994464819411, 51.5067557678594, 51.4977275537675, 51.4988718377984, 51.4974782970503, 51.5137643501102, 51.5204498816501, 51.5213788858189, 51.5144343479237, 51.5212383455566, 51.5088225370868, 51.512547894056, 51.5144758355252, 51.5218865924773, 51.5218808497196, 51.5152330574081, 51.4836680563637, 51.4885211991595, 51.486842412489, 51.48845301363455)
  } else {
  	if (is.null(longitude) | is.null(latitude))
  		stop("Both longitude and latitude arguments must be used, or alternatively both arguments must be set to NULL to use default values")
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
#' @param longitude the longitudinal positions of the potential sources.
#' @param latitude the latitudinal positions of the potential sources.
#'
#' @export
#' @examples
#' geoDataSource()

geoDataSource <- function(source_longitude=NULL, source_latitude=NULL) {
    
  # generate dummy data if none read in
  if (is.null(source_longitude) & is.null(source_latitude)) {
      source_longitude <- c(-0.1, -0.14, -0.105, -0.08, -0.08)
      source_latitude <- c(51.495, 51.5, 51.52, 51.515, 51.49)
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
#' This function can be used to generate parameters in the format required by other Rgeoprofile functions. Any parameter value can be specified as an input argument to this function. Alternatively, if data is input as an argument then some parameters can take default values directly from the data.
#'
#' @param data observations in the format defined by geoData().
#' @param sigma_mean the mean of the prior on sigma (sigma = standard deviation of the dispersal distribution).
#' @param sigma_var the variance of the prior on sigma.
#' @param sigma_squared_shape as an alternative to defining the prior mean and variance of sigma, it is possible to directly define the parameters of the inverse-gamma prior on sigma^2. This is the shape parameter of the inverse-gamma prior.
#' @param sigma_squared_rate this is the rate parameter of the inverse-gamma prior on sigma^2.
#' @param priorMean_longitude the position of the prior mean on source locations in degrees longitude.
#' @param priorMean_latitude the position of the prior mean on source locations in degrees latitude.
#' @param tau the standard deviation of the normal prior on source locations, i.e. how far we expect sources to lie from the centre. Leave as NULL to use default value, in which case tau is set equal to the maximum distance of any observation from the prior mean.
#' @param alpha_shape shape parameter of the gamma prior on the parameter alpha.
#' @param alpha_rate rate parameter of the gamma prior on the parameter alpha.
#' @param chains number of MCMC chains to use in the burn-in step.
#' @param burnin number of burn-in iterations to be discarded at start of MCMC.
#' @param samples number of sampling iterations. These iterations are used to generate final posterior distribution.
#' @param burnin_printConsole how frequently (in iterations) to report progress to the console during the burn-in phase.
#' @param samples_printConsole how frequently (in iterations) to report progress to the console during the sampling phase.
#' @param longitude_minMax vector containing minimum and maximum longitude over which to generate geoprofile.
#' @param latitude_minMax vector containing minimum and maximum latitude over which to generate geoprofile.
#' @param longitude_cells number of cells in the final geoprofile (longitude direction). Higher values generate smoother distributions, but take longer to run.
#' @param latitude_cells number of cells in the final geoprofile (latitude direction). Higher values generate smoother distributions, but take longer to run.
#'
#' @export
#' @examples
#' myData <- geoData()
#' geoParams(myData, sigma_var=1)

geoParams <- function(data=NULL, sigma_mean=1, sigma_var=NULL, sigma_squared_shape=NULL, sigma_squared_rate=NULL, priorMean_longitude=NULL, priorMean_latitude=NULL, tau=NULL, alpha_shape=0.1, alpha_rate=0.1, chains=10, burnin=500, samples=5000, burnin_printConsole=100, samples_printConsole=1000, longitude_minMax=NULL, latitude_minMax=NULL, longitude_cells=500, latitude_cells=500, guardRail=0.05) {
    
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
			longitude_minMax <- priorMean_longitude + c(-0.1,0.1)
		if (is.null(latitude_minMax))
			latitude_minMax <- priorMean_latitude + c(-0.1,0.1)
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
#' @param silent Wwether to report passing check to console.
#'
#' @export
#' @examples
#' myData <- geoData()
#' geoDataCheck(myData)

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
#' myData <- geoData()
#' myParams <- geoParams(myData, sigma_var=1)
#' geoParamsCheck(myParams)

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
#' @param sigma a vector of posterior draws of sigma. Leave as NULL to plot prior only.
#' @param plotMax maximum x-axis range to plot. Leave as NULL to use default settings.
#'
#' @export
#' @examples
#' myData <- geoData()
#' myParams <- geoParams(myData, sigma_var=1)
#' geoPlotSigma(myParams)

geoPlotSigma <- function(params, sigma=NULL, plotMax=NULL) {
  
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
  
  # default plotMax based on extent of prior distribution AND the extent of posterior draws if available
  if (is.null(plotMax)) {
    plotMax <- sigma_mean+3*sqrt(sigma_var)
    if (!is.null(sigma)) {
      plotMax <- max(plotMax, 2*max(sigma,na.rm=TRUE))
    }
  }
  
  # produce prior distribution
  sigma_vec <- seq(0,plotMax,l=501)
  sigma_prior <- dRIG(sigma_vec,alpha,beta)
  
  # plot prior and overlay density of posterior draws if used
  if (is.null(sigma)) {
    plot(sigma_vec, sigma_prior, type='l', xlab='sigma (km)', ylab='probability density', main='')
    legend(x='topright', legend='prior', lty=1)
  } else {
    sigma_posterior <- density(sigma,from=0,to=plotMax)
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
#' MCMC under Rgeoprofile model
#'
#' This function carries out the main MCMC under the Rgeoprofile model.
#'
#' @param data input data in the format defined by geoData().
#' @param params input parameters in the format defined by geoParams().
#'
#' @export
#' @examples
#' myData <- geoData()
#' myParams <- geoParams(myData, sigma_var=1)
#' myMCMC <- geoMCMC(myData, myParams)
#' geoPlotSigma(myParams, myMCMC$sigma)

geoMCMC <- function(data, params) {

  # check that data and parameters in correct format
  geoDataCheck(data)
  geoParamsCheck(params)
  cat("\n")
  
  # transform data and map limits to cartesian coordinates relative to centre of prior. After transformation data are defined relative to point 0,0 (i.e. the origin represents the centre of the prior)
  data_cartesian <-latlon_to_cartesian(params$model$priorMean_latitude, params$model$priorMean_longitude, data$latitude, data$longitude)
  limits_cartesian <-latlon_to_cartesian(params$model$priorMean_latitude, params$model$priorMean_longitude, params$output$latitude_minMax, params$output$longitude_minMax)
  
  # add these cartesian coordinates to data and params objects before feeding into C++ function
  data$x <- data_cartesian$x
  data$y <- data_cartesian$y
  params$output$x_minMax <- limits_cartesian$x
  params$output$y_minMax <- limits_cartesian$y
  
  # if using fixed sigma model then change alpha and beta from NULL to -1. This value will be ignored, but needs to be numeric before feeding into the C++ function.
  if (params$model$sigma_var==0) {
  	params$model$sigma_squared_shape <- -1
  	params$model$sigma_squared_rate <- -1
  }
    
  # carry out MCMC using efficient C++ function
  rawOutput <- C_geoMCMC(data, params)
  
  # extract raw draws and check that at least one posterior draw in chosen region
  surface_raw <- matrix(unlist(rawOutput$geoSurface), params$output$latitude_cells, byrow=TRUE)
  if (all(surface_raw==0))
      stop('chosen lat/long window contains no posterior draws')
  
  # get size of each cell
  cells_x <- params$output$longitude_cells
  cells_y <- params$output$latitude_cells
  cellSize_x <- diff(limits_cartesian$x)/cells_x
  cellSize_y <- diff(limits_cartesian$y)/cells_y
  
  # temporarily add guard rail to surface to avoid Fourier series bleeding round edges
  railSize_x <- cells_x
  railSize_y <- cells_y
  railMat_x <- matrix(0,cells_y,railSize_x)
  railMat_y <- matrix(0,railSize_y,cells_x+2*railSize_x)
  
  surface_normalised <- surface_raw/sum(surface_raw)
  surface_normalised <- cbind(railMat_x, surface_normalised, railMat_x)
  surface_normalised <- rbind(railMat_y, surface_normalised, railMat_y)
  
  # calculate Fourier transform of posterior surface
  f1 = fftw2d(surface_normalised)
  
  # produce surface that kernel will be calculated over
  kernel_x <- cellSize_x * c(0:floor(ncol(surface_normalised)/2), floor((ncol(surface_normalised)-1)/2):1)
  kernel_y <- cellSize_y * c(0:floor(nrow(surface_normalised)/2), floor((nrow(surface_normalised)-1)/2):1)
  kernel_x_mat <- outer(rep(1,length(kernel_y)), kernel_x)
  kernel_y_mat <- outer(kernel_y, rep(1,length(kernel_x)))
  kernel_s_mat <- sqrt(kernel_x_mat^2+kernel_y_mat^2)
  
  # set lambda (bandwidth) increment size based on cell size
  lambda_step <- min(cellSize_x, cellSize_y)/5
  
  # loop through range of values of lambda
  cat('Smoothing posterior surface\n')
  flush.console()
  logLike <- -Inf
  for (i in 1:100) {
      
    # calculate Fourier transform of kernel
    lambda <- lambda_step*i
    kernel <- dts(kernel_s_mat,df=3,scale=lambda)
    f2 = fftw2d(kernel)
    
    # combine Fourier transformed surfaces and take inverse. f4 will ultimately become the main surface of interest.
    f3 = f1*f2
    f4 = Re(fftw2d(f3,inverse=T))/length(surface_normalised)
    
    # subtract from f4 the probability density of each point measured from itself. In other words, move towards a leave-one-out kernel density method
    f5 <- f4 - surface_normalised*dts(0,df=3,scale=lambda)
    f5[f5<0] <- 0
    f5 <- f5/sum(f4)
    
    # calculate leave-one-out log-likelihood at each point on surface
    f6 <- surface_normalised*log(f5)
    
    # break if total log-likelihood is at a maximum
    if (sum(f6,na.rm=T)<logLike)
      break()
    logLike <- sum(f6,na.rm=T)

  }
  cat(paste('maximum likelihood lambda = ',round(lambda,3),sep=''))
  
  # remove guard rail
  f4 <- f4[,(railSize_x+1):(ncol(f4)-railSize_x)]
  f4 <- f4[(railSize_y+1):(nrow(f4)-railSize_y),]
  
  # produce prior matrix. Note that each cell of this matrix contains the probability density at that point multiplied by the size of that cell, meaning the total sum of the matrix from -infinity to +infinity would equal 1. As the matrix is limited to the region specified by the limits, in reality this matrix will sum to some value less than 1.
  x_mids <- seq(limits_cartesian$x[1], limits_cartesian$x[2], l=ncol(f4)+1)[-1] - cellSize_x/2
  y_mids <- seq(limits_cartesian$y[1], limits_cartesian$y[2], l=nrow(f4)+1)[-1] - cellSize_y/2
  x_mids_mat <- outer(rep(1,cells_y),x_mids)
  y_mids_mat <- outer(y_mids,rep(1,cells_x))
  
  priorMat <- dnorm(x_mids_mat,sd=params$model$tau)*dnorm(y_mids_mat,sd=params$model$tau)*(cellSize_x*cellSize_y)
  
  # finalise output format
  output <- list()
  
  # sigma
  output$sigma <- rawOutput$sigma
  
  # alpha
  alpha <- rawOutput$alpha
  output$alpha <- alpha
  
  # combine prior surface with stored posterior surface (the prior never fully goes away under a DPM model)
  output$surface_raw <- surface_raw
  n <- length(data$longitude)
  output$surface <-  f4 + priorMat*mean(alpha/(alpha+n))
  
  # posterior allocation
  allocation <- matrix(unlist(rawOutput$allocation),n,byrow=T)
  allocation <- data.frame(allocation/params$MCMC$samples)
  names(allocation) <- paste("group",1:ncol(allocation),sep="")
  output$allocation <- allocation
  
  return(output)
}

#------------------------------------------------
#' kernel density smoothing of posterior distribution
#'
#' Can be used to perform kernel density smoothing of the "surface_raw" object output from geoMCMC(). Note that geoMCMC() performs this smoothing already using a maximum-likelihood estimate of lambda and outputs it to the "surface" object, and so this function is only needed when a custom level of smoothing is required. The kernel used is a Student's t distribution with a user-defined scale (1 = ordinary Student's t distribution) and degrees of freedom.
#'
#' @param data input data in the format defined by geoData().
#' @param params input parameters in the format defined by geoParams().
#' @param MCMCoutput stored output of the MCMC obtained by running geoMCMC().
#' @param lambda scale of smoothing kernel, relative to an ordinary Student's t distribution in which lambda=1. When appled to geoMCMC() output this parameter is in units of km.
#' @param df degrees of freedom of Student's t smoothing kernel.
#'
#' @export
#' @examples
#' myData <- geoData()
#' myParams <- geoParams(myData, sigma_var=1)
#' myMCMC <- geoMCMC(myData, myParams)
#' mySurface <- geoSmooth(myData, myParams, myMCMC, lambda=0.5, df=2)

geoSmooth <- function(data, params, MCMCoutput, lambda=1, df=3) {

  # get size of each cell
  limits_cartesian <-latlon_to_cartesian(params$model$priorMean_latitude, params$model$priorMean_longitude, params$output$latitude_minMax, params$output$longitude_minMax)
  cells_x <- params$output$longitude_cells
  cells_y <- params$output$latitude_cells
  cellSize_x <- diff(limits_cartesian$x)/cells_x
  cellSize_y <- diff(limits_cartesian$x)/cells_y

  # temporarily add guard rail to surface to avoid Fourier series bleeding round edges
  railSize_x <- cells_x
  railSize_y <- cells_y
  railMat_x <- matrix(0,cells_y,railSize_x)
  railMat_y <- matrix(0,railSize_y,cells_x+2*railSize_x)
  
  surface_normalised <- MCMCoutput$surface_raw/sum(MCMCoutput$surface_raw, na.rm=TRUE)
  surface_normalised <- cbind(railMat_x, surface_normalised, railMat_x)
  surface_normalised <- rbind(railMat_y, surface_normalised, railMat_y)
  
  # calculate Fourier transform of posterior surface
  f1 = fftw2d(surface_normalised)
  
  # produce surface that kernel will be calculated over
  kernel_x <- cellSize_x * c(0:floor(ncol(surface_normalised)/2), floor((ncol(surface_normalised)-1)/2):1)
  kernel_y <- cellSize_y * c(0:floor(nrow(surface_normalised)/2), floor((nrow(surface_normalised)-1)/2):1)
  kernel_x_mat <- outer(rep(1,length(kernel_y)), kernel_x)
  kernel_y_mat <- outer(kernel_y, rep(1,length(kernel_x)))
  kernel_s_mat <- sqrt(kernel_x_mat^2+kernel_y_mat^2)
  
  # calculate Fourier transform of kernel
  kernel <- dts(kernel_s_mat,df=3,scale=lambda)
  f2 = fftw2d(kernel)
    
  # combine Fourier transformed surfaces and take inverse. f4 will ultimately become the main surface of interest.
  f3 = f1*f2
  f4 = Re(fftw2d(f3,inverse=T))/length(surface_normalised)
    
  # remove guard rail
  f4 <- f4[,(railSize_x+1):(ncol(f4)-railSize_x)]
  f4 <- f4[(railSize_y+1):(nrow(f4)-railSize_y),]
  
  # produce prior matrix. Note that each cell of this matrix contains the probability density at that point multiplied by the size of that cell, meaning the total sum of the matrix from -infinity to +infinity would equal 1. As the matrix is limited to the region specified by the limits, in reality this matrix will sum to some value less than 1.
  x_mids <- seq(limits_cartesian$x[1], limits_cartesian$x[2], l=ncol(f4)+1)[-1] - cellSize_x/2
  y_mids <- seq(limits_cartesian$y[1], limits_cartesian$y[2], l=nrow(f4)+1)[-1] - cellSize_y/2
  x_mids_mat <- outer(rep(1,cells_y),x_mids)
  y_mids_mat <- outer(y_mids,rep(1,cells_x))
  
  priorMat <- dnorm(x_mids_mat,sd=params$model$tau)*dnorm(y_mids_mat,sd=params$model$tau)*(cellSize_x*cellSize_y)
  
  # combine prior surface with stored posterior surface (the prior never fully goes away under a DPM model)
  n <- length(data$longitude)
  alpha <- MCMCoutput$alpha
  output <-  f4 + priorMat*mean(alpha/(alpha+n))
  
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
#' myData <- geoData()
#' myParams <- geoParams(myData, sigma_var=1)
#' myMCMC <- geoMCMC(myData, myParams)
#' myMCMC$profile <- geoProfile(myMCMC$surface)

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
#' @param MCMCoutput output generated from an MCMC.
#'
#' @export
#' @examples
#' geoPlotAllocation(MCMCoutput)

geoPlotAllocation <- function(allocation, colours="default", barBorderCol="white", barBorderWidth=0.25, mainBorderCol="black", mainBorderWidth=2, yTicks_on=TRUE, yTicks=seq(0,1,0.2), xlab="", ylab="posterior allocation", mainTitle="", names=NA, names_size=1, xTicks_on=FALSE, xTicks_size=1, orderBy="group") {
    
    # check that orderBy is either 'group' or 'probability'
    if (!(orderBy%in%c("group","probability")))
        stop("orderBy must equal 'group' or 'probability'")
    
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
#' Uses ggmap to download a map for the region specified by params. Various elements can be overlaid onto this map, including raw observation data, potential source locations, and a final geoprofile. Plotting options such as the type of map, the size and colour of points, and the contours of the geoprofile can be controlled.
#'
#' @param params parameters list in the format defined by geoParams().
#' @param data data object in the format defined by geoData().
#' @param source potential sources object in the format defined by geoDataSource().
#' @param surface a surface to overlay onto the map, typically a geoprofile produced from the output of geoMCMC().
#' @param zoom some text.
#' @param mapSource which online source to use when downloading the map. Options include Google Maps ("google"), OpenStreetMap ("osm"), Stamen Maps ("stamen") and CloudMade maps ("cloudmade").
#' @param params mapType the specific type of map of map to plot. Options available are "terrain", "terrain-background", "satellite", "roadmap" and "hybrid" (google maps), "terrain", "watercolor" and "toner" (stamen maps) or a positive integer for cloudmade maps (see ?get_cloudmademap from the package ggmap).
#' @param transparency some text.
#' @param plotContours some text.
#' @param breakPercent some text.
#' @param contourCols some text.
#' @param crimeCex some text.
#' @param crimeCol some text.
#' @param crimeBorderCol some text.
#' @param crimeBorderWidth some text.
#' @param sourceCex some text.
#' @param sourceCol some text.
#'
#' @export
#' @examples
#' geoQuickPlot(surface)

geoPlotMap <- function(params, data=NULL, source=NULL, surface=NULL, zoom="auto", mapSource="google", mapType="hybrid", transparency=0.6, plotContours=TRUE, breakPercent=seq(0,100,l=11), contourCols= c("red","orange","yellow","white"), crimeCex=1.5, crimeCol='red', crimeBorderCol='white', crimeBorderWidth=0.5, sourceCex=1.5, sourceCol='blue') {
    
    # check that inputs make sense
    geoParamsCheck(params)
    if (!is.null(data))
	    geoDataCheck(data)
    
    # if zoom=="auto" then set zoom level based on params
    if (zoom=="auto")
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
		myMap <- myMap + geom_tile(aes(x=x,y=y,fill=cut), alpha=transparency, data=df_noNA)
		myMap <- myMap + scale_fill_manual(name="Hitscore\npercentage", values=rev(geoCols(nbcol)))

		# add contours
		if (plotContours) {
			myMap <- myMap + stat_contour(aes(x=x,y=y,z=z), colour="grey50", breaks=breakPercent/100*length(surface), size=0.3, alpha=transparency, data=df)
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
#' Calculate hitscores of the potential sources based on the final geoprofile surface.
#'
#' @param params some text
#' @param source_data some text
#' @param surface some text
#'
#' @export
#' @examples
#' geoReportHitscores(params,source_data,surface)

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
		hitscores = hitscoremat2[msourcey,msourcex]
	}
	hit_output <<- cbind(sources,hitscores)
	print(hit_output)

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
#' geoPlotLorenz()

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
#' @param crime1 numerical index of first crime
#' @param crime2 numerical index of second crime
#' @param coallocation_matrix matrix of coallocations between all observations, as produced the "allocation" output of the function geoMCMC()
#' @param offset vertical offset of second line to ensure readability
#' @param plot.graph whether to plot the graph
#'
#' @export
#' @examples
#' prob_coallocation()

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
#' Calculate and plot ring search strategy
#'
#' Calculates hit scores for a ring-search strategy (ie searching in an expanding radius out from the crimes).
#' Also plots the crimes and sources with merged polygons showing these (merged and clipped) rings
#'
#' @param data some text
#' @param source some text
#' @param buffer_radii some text
#'
#' @export
#' @examples
#' ringHS()

ringHS <- function(params,crime_data, source_data, buffer_radii=c(1000,2000,5000))
	{
        library(RgoogleMaps)
        library(rgeos)
        # function for calculating UTM zone from mean of longitude of crimes
		long2UTM <- function(long)
			{
				(floor((long + 180)/6) %% 60) + 1
			}

		# calculate UTM zone
		my_UTM <- long2UTM(mean(crime_data$longitude))
		
		# set projections for ease of use
		my_crs_long_lat <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
		my_crs_utm <- paste("+proj=utm +zone=",my_UTM," ellps=WGS84",sep="")

		##load crime and source data and convert to spatial object with projection of choice (UTM and lonlat)
		crimes <- cbind(crime_data$longitude, crime_data$latitude)
		crimes_lonlat <- SpatialPointsDataFrame(coords = as.matrix(crimes), data = as.data.frame(crimes),proj4string = CRS(my_crs_long_lat))
		sources <- cbind(source_data$source_longitude, source_data $source_latitude)
		sources_lonlat <- SpatialPointsDataFrame(coords = as.matrix(sources), data = as.data.frame(sources),proj4string = CRS(my_crs_long_lat))
		# plot as check
		# plot(crimes_lonlat)
		# plot(sources_lonlat)
		
		# use params to calculate area of map
		lonMin <- params$output$longitude_minMax[1]
		lonMax <- params$output$longitude_minMax[2]
		latMin <- params$output$latitude_minMax[1]
		latMax <- params$output$latitude_minMax[2]
		bounds=matrix(c(lonMin, latMin,
				lonMin, latMax,
				lonMax, latMax,
				lonMax,  latMin,
				lonMin, latMin
				),
				ncol=2, byrow=TRUE)
		# plot as check		
		# plot(bounds,type="l")
		# create polygon as spatial object of bounding box as both longlat and UTM
		bounds_lonlat <- SpatialPointsDataFrame(coords = as.matrix(bounds), data = as.data.frame(bounds),proj4string = CRS(my_crs_long_lat))
		bounds_UTM <- spTransform(bounds_lonlat, CRS(my_crs_utm))
        
		# turn bounding box into a polygon, and then a spatial polygon with UTM projection to match merged circles
		p1=Polygon(bounds_UTM)
		bounds_polygon_utm = SpatialPolygons(list(Polygons(list(p1), ID = "a")), proj4string=CRS(my_crs_utm))
		# plot(bounds_polygon_utm, axes = TRUE)

		# area in metres square of map (for hit score calculations later)
		map_area_m_sq <- gArea(bounds_polygon_utm)
		
		#convert to UTM so we can use m for gBuffer function, width is radius in m
		crimes_UTM <- spTransform(crimes_lonlat, CRS(my_crs_utm))
		# plot(crimes_UTM)
		sources_UTM <- spTransform(sources_lonlat, CRS(my_crs_utm))
		# plot(sources_UTM)

		# extract number of sources and calculate for each the distance to the nearest crime (in m)
		no_sources <- dim(sources_lonlat)[1]
		min_dists <- rep(NA,no_sources)
		for(i in 1:no_sources)
			{
				min_dists[i] <- min(spDistsN1(crimes_lonlat, sources_lonlat[i,1],longlat=TRUE))*1000
			}
		# min_dists

		# calculate ring hit scores for each source
		stored_results <- matrix(rep(NA, no_sources*4),ncol=4)
		stored_buffers <- list()
		colnames(stored_results) <- c("ring_lon","ring_lat","merged_buffer_area_m2","ring_hs")
		for(source_number_to_check in 1: no_sources)
			{
				# create buffer
				b <- gBuffer(crimes_UTM, byid = FALSE, width = min_dists[source_number_to_check])
				# plot(b)
				#merge the circles
				b2<- gUnaryUnion(b)
				clip <-gIntersection(b2, bounds_polygon_utm, byid=TRUE, drop_lower_td=TRUE)
				
				#calculate the area
				merged_area <- gArea(clip)
				ring_hs_for_this_source <- merged_area/map_area_m_sq
				stored_results[source_number_to_check,] <- c(sources[source_number_to_check,],merged_area, ring_hs_for_this_source)
				stored_buffers[source_number_to_check] <- clip
			}
		
		# plot
		# create list to store contours for plotting on map
		stored_contours <- list()
		for(contour_number in 1:length(buffer_radii))
			{
				# create buffer
				b <- gBuffer(crimes_UTM, byid = FALSE, width = buffer_radii[contour_number])
				# plot(b)
				#merge the circles
				b2<- gUnaryUnion(b)
				clip <-gIntersection(b2, bounds_polygon_utm, byid=TRUE, drop_lower_td=TRUE)
				
				clip_lonlat <- spTransform(clip,CRS(my_crs_long_lat))

				stored_contours[contour_number] <- clip_lonlat
			}

		# download map
		MyMap <- GetMap.bbox(params$output$longitude_minMax, params$output$latitude_minMax,maptype="roadmap")
		MyMap <- GetMap.bbox(params$output$longitude_minMax, params$output$latitude_minMax,maptype="roadmap",destfile="MyTile.png")

		# plot map
		quartz("ring search")
		PlotOnStaticMap(MyMap)
		
		gp.colors <- colorRampPalette(my_reds)

		ringCols <-  gp.colors(length(stored_contours))
		ringColsTransp = AddAlpha(ringCols,0.05)
		
		# plot contours
		for(cc in length(stored_contours):1)
			{
				PlotPolysOnStaticMap(MyMap, stored_contours[[cc]],add=TRUE,col=ringColsTransp[cc],border= "black",lwd=1)
			}
		# plot bounding box
		box_to_plot <- spTransform(bounds_polygon_utm,CRS(my_crs_long_lat))
		PlotPolysOnStaticMap(MyMap, box_to_plot,add=TRUE,col=NULL,border="black",lwd=2)
		# add crimes and sources
		PlotOnStaticMap(MyMap,crime_data$latitude,crime_data$longitude, pch=16,add=TRUE)
		PlotOnStaticMap(MyMap,source_data$source_latitude,source_data$source_longitude, pch=15,col="red",add=TRUE)
	
		ring_table<-data.frame(stored_results)
		ring_hs <- ring_table[,4]
		ring_areas <- ring_table[,3]
		
		# return results
		# with polygons
		# ring_output <- list(ring_table=ring_table,ring_areas=ring_areas,ring_hs=ring_hs,my_UTM= my_UTM,map_area_m_sq=map_area_m_sq,stored_buffers=stored_buffers)
		# without polygons
		ring_output <- list(ring_table=ring_table,ring_areas=ring_areas,ring_hs=ring_hs,my_UTM=my_UTM,map_area_m_sq=map_area_m_sq)
		return(ring_output)
}
		
#------------------------------------------------
