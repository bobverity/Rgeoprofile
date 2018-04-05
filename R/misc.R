
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
#' Convert lat-lon to bearing
#'
#' Calculate bearing and great circle distance between an origin and one or more destination points
#'
#' @param origin_lat latitude of origin point
#' @param origin_lon longitude of origin point
#' @param dest_lat latitude of destination point
#' @param dest_lon longitude of destination point

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

dts <- function(x, df, scale=1, log=FALSE) {
  ret <- lgamma((df+1)/2)-lgamma(df/2)-0.5*log(pi*df*scale^2) - ((df+1)/2)*log(1 + x^2/(df*scale^2))
  if (!log) { ret <- exp(ret) }
  return(ret)
}

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
# Draw from normal distribution converted to spherical coordinate system. Points are first drawn from an ordinary cartesian 2D normal distribution. The distances to points are then assumed to be great circle distances, and are combined with a random bearing from the point {centre_lat, centre_lon} to produce a final set of lat/lon points. Note that this is not a truly spherical normal distribution, as the domain of the distribution is not the sphere - rather it is a transformation from one coordinate system to another that is satisfactory when the curvature of the sphere is not severe.
# (not exported)

rnorm_sphere <- function(n, centre_lat, centre_lon, sigma) {
  x <- rnorm(n,sd=sigma)
  y <- rnorm(n,sd=sigma)
  output <- cartesian_to_latlon(centre_lat, centre_lon, x, y)
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
#' Produce a smooth surface using 2D kernel density smoothing
#'
#' Takes lon/lat coordinates, bins in two dimensions, and smooths using kernel density smoothing. Kernel densities are computed using the fast Fourier transform method, which is many times faster than simple summation when using a large number of points. Each Kernel is student's-t distributed with 3 degrees of freedom, and scaled by the bandwidth lambda. If lambda is set to \code{NULL} then the optimal value of lambda is chosen automatically using the leave-one-out maximum likelihood method.
#'
#' @param longitude longitude of input points
#' @param latitude latitude of input points
#' @param breaks_lon positions of longitude breaks
#' @param breaks_lat positions of latitude breaks
#' @param lambda bandwidth to use in posterior smoothing. If NULL then optimal bandwidth is chosen automatically by maximum-likelihood.
#'
#' @export
#' @examples
#' # TODO

geoSmooth <- function(longitude, latitude, breaks_lon, breaks_lat, lambda=NULL) {
  
  # get properties of cells in each dimension
  cells_lon <- length(breaks_lon) - 1
  cells_lat <- length(breaks_lat) - 1
  centre_lon <- mean(breaks_lon)
  centre_lat <- mean(breaks_lat)
  cellSize_lon <- diff(breaks_lon[1:2])
  cellSize_lat <- diff(breaks_lat[1:2])
  
  # bin lon/lat values in two dimensions and check that at least one value in chosen region
  surface_raw <- bin2D(longitude, latitude, breaks_lon, breaks_lat)$z
  if (all(surface_raw==0)) { stop('chosen lat/long window contains no posterior draws') }
  
  # temporarily add guard rail to surface to avoid Fourier series bleeding round edges
  railSize_lon <- cells_lon
  railSize_lat <- cells_lat
  railMat_lon <- matrix(0, cells_lat, railSize_lon)
  railMat_lat <- matrix(0, railSize_lat, cells_lon + 2*railSize_lon)
  
  surface_normalised <- surface_raw/sum(surface_raw)
  surface_normalised <- cbind(railMat_lon, surface_normalised, railMat_lon)
  surface_normalised <- rbind(railMat_lat, surface_normalised, railMat_lat)
  
  # calculate Fourier transform of posterior surface
  f1 = fftw2d(surface_normalised)
  
  # calculate x and y size of one cell in cartesian space. Because of transformation, this size will technically be different for each cell, but use centre of space to get a middling value
  cellSize_trans <- latlon_to_cartesian(centre_lat, centre_lon, centre_lat + cellSize_lat, centre_lon + cellSize_lon)
  cellSize_trans_lon <- cellSize_trans$x
  cellSize_trans_lat <- cellSize_trans$y
  
  # produce surface over which kernel will be calculated. This surface wraps around in both x and y (i.e. the kernel is actually defined over a torus).
  kernel_lon <- cellSize_trans_lon * c(0:floor(ncol(surface_normalised)/2), floor((ncol(surface_normalised) - 1)/2):1)
  kernel_lat <- cellSize_trans_lat * c(0:floor(nrow(surface_normalised)/2), floor((nrow(surface_normalised) - 1)/2):1)
  kernel_lon_mat <- outer(rep(1,length(kernel_lat)), kernel_lon)
  kernel_lat_mat <- outer(kernel_lat, rep(1,length(kernel_lon)))
  kernel_s_mat <- sqrt(kernel_lon_mat^2 + kernel_lat_mat^2)
  
  # set lambda (bandwidth) range to be explored
  if (is.null(lambda)) {
    lambda_step <- min(cellSize_trans_lon, cellSize_trans_lat)/5
    lambda_vec <- lambda_step*(1:100)
  } else {
    lambda_vec <- lambda
  }
  
  # loop through range of values of lambda
  cat('Smoothing posterior surface')
  flush.console()
  logLike <- -Inf
  for (i in 1:length(lambda_vec)) {
    
    # print dots to screen
    if (i>1) {
      cat(".")
      flush.console()
    }
    
    # calculate Fourier transform of kernel
    lambda_this <- lambda_vec[i]
    kernel <- dts(kernel_s_mat, df=3, scale=lambda_this)
    f2 = fftw2d(kernel)
    
    # combine Fourier transformed surfaces and take inverse. f4 will ultimately become the main surface of interest.
    f3 = f1*f2
    f4 = Re(fftw2d(f3,inverse=T))/length(surface_normalised)
    
    # subtract from f4 the probability density of each point measured from itself. In other words, move towards a leave-one-out kernel density method
    f5 <- f4 - surface_normalised*dts(0, df=3, scale=lambda_this)
    f5[f5<0] <- 0
    f5 <- f5/sum(f4)
    
    # calculate leave-one-out log-likelihood at each point on surface
    f6 <- surface_normalised*log(f5)
    
    # break if total log-likelihood is at a local maximum
    if (sum(f6,na.rm=T)<logLike) { break() }
    
    # otherwise update logLike
    logLike <- sum(f6,na.rm=T)
  }
  
  # report chosen value of lambda
  if (is.null(lambda)) {
    cat(paste('\nmaximum likelihood lambda = ', round(lambda_this,3), sep=''))
  }
  
  # remove guard rail
  f4 <- f4[,(railSize_lon+1):(ncol(f4)-railSize_lon)]
  f4 <- f4[(railSize_lat+1):(nrow(f4)-railSize_lat),]
  
  # return surface
  return(f4)
}
