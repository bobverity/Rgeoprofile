
#------------------------------------------------
#' @title Import file
#'
#' @description Import file from the inst/extdata folder of this package
#' 
#' @param name name of file
#'
#' @export

rgeoprofile_file <- function(name) {
  
  # load file from inst/extdata folder
  name_full <- system.file("extdata/", name, package = 'RgeoProfile', mustWork = TRUE)
  ret <- readRDS(name_full)
  
  # return
  return(ret)
}

#------------------------------------------------
# replace NULL value with default
#' @noRd
define_default <- function(x, default_value) {
  if (is.null(x)) {
    x <- default_value
  }
  return(x)
}

#------------------------------------------------
# simple zero-padding function. Not robust to e.g. negative numbers
#' @noRd
zero_pad_simple <- function(x, n = 3) {
  ret <- mapply(function(x) {
                  paste0(paste0(rep(0,n-nchar(x)), collapse = ""), x, collapse = "")
                }, x)
  return(ret)
}

# -----------------------------------
# ask user a yes/no question. Return TRUE/FALSE
#' @noRd
user_yes_no <- function(x="continue? (Y/N): ") {
  userChoice <- NA
  while (!userChoice %in% c("Y", "y" ,"N", "n")) {
    userChoice <- readline(x)
  }
  return(userChoice %in% c("Y", "y"))
}

# -----------------------------------
# draw from Dirichlet distribution
#' @noRd
rdirichlet <- function (alpha_vec) {
  Y <- rgamma(length(alpha_vec), shape = alpha_vec, scale = 1)
  output <- Y/sum(Y)
  return(output)
}

# -----------------------------------
# takes matrix as input, converts to list format for use within Rcpp code
#' @noRd
mat_to_rcpp <- function(x) {
  return(split(x, f=1:nrow(x)))
}

# -----------------------------------
# takes list format returned from Rcpp and converts to matrix
#' @noRd
rcpp_to_mat <- function(x) {
  ret <- matrix(unlist(x), nrow=length(x), byrow=TRUE)
  return(ret)
}

#------------------------------------------------
# calls C++ implementation of the Hungarian algorithm for finding best matching
# in a linear sum assigment problem. This is function is used in testing.
#' @noRd
call_hungarian <- function(x) {
  args <- list(cost_mat = mat_to_rcpp(x))
  call_hungarian_cpp(args)
}

#------------------------------------------------
# return 95% quantile
#' @noRd
quantile_95 <- function(x) {
  ret <- quantile(x, probs=c(0.025, 0.5, 0.975))
  names(ret) <- c("Q2.5", "Q50", "Q97.5")
  return(ret)
}

#------------------------------------------------
# sum logged values without underflow, i.e. do log(sum(exp(x)))
#' @noRd
log_sum <- function(x) {
  if (all(is.na(x))) {
    return(rep(NA, length(x)))
  }
  x_max <- max(x, na.rm = TRUE)
  ret <- x_max + log(sum(exp(x-x_max)))
  return(ret)
}

#------------------------------------------------
# return p-value of Geweke's diagnostic convergence statistic, estimated from package coda
#' @noRd
geweke_pvalue <- function(x) {
  ret <- 2*pnorm(abs(geweke.diag(x)$z), lower.tail=FALSE)
  return(ret)
}

#------------------------------------------------
# check that geweke p-value non-significant on values x[1:n]
#' @noRd
test_convergence <- function(x, n) {
  if (n==1) {
    return(FALSE)
  }
  g <- geweke_pvalue(mcmc(x[1:n]))
  ret <- (g>0.01)
  if (is.na(ret)) {
    ret <- TRUE;
  }
  return(ret)
}

#------------------------------------------------
# update progress bar
# (not exported)
#' @noRd
update_progress <- function(pb_list, name, i, max_i) {
  setTxtProgressBar(pb_list[[name]], i)
  if (i == max_i) {
    close(pb_list[[name]])
  }
}

#------------------------------------------------
#' @title Draw from spherical normal distribution
#'
#' @description Draw from normal distribution converted to spherical coordinate
#'   system. Points are first drawn from an ordinary cartesian 2D normal
#'   distribution. The distances to points are then assumed to be great circle
#'   distances, and are combined with a random bearing from the point
#'   {centre_lat, centre_lon} to produce a final set of lat/lon points.
#'
#' @param n The number of points to draw
#' @param centre_lon The mean longitude of the normal distribution
#' @param centre_lat The mean latitude of the normal distribution
#' @param sigma The standard deviation of the normal distribution
#'
#' @export
#' @examples
#' rnorm_sphere(n = 100, centre_lat = 0, centre_lon = 0, sigma = 1)

rnorm_sphere <- function(n, centre_lon, centre_lat, sigma = 1) {
  
  # draw points centred at zero
  x <- rnorm(n, sd = sigma)
  y <- rnorm(n, sd = sigma)
  
  # calculate angle and euclidian distance of all points from origin. Angles are
  # in degrees relative to due north
  d <- sqrt(x^2 + y^2)
  theta <- atan2(x, y)*360/(2*pi)
  
  # get lon/lat relative to origin
  ret <- bearing_to_lonlat(centre_lon, centre_lat, theta, d)
  
  return(ret)
}

#------------------------------------------------
#' @title Get spatial coordinate given an origin, a great circle distance and a
#'   bearing
#'
#' @description Calculate destination lat/lon given an origin, a great circle
#'   distance of travel, and a bearing.
#'
#' @param origin_lon The origin longitude
#' @param origin_lat The origin latitude
#' @param bearing The angle in degrees relative to due north
#' @param gc_dist The great circle distance in (km)
#'
#' @export
#' @examples
#' # one degree longitude is approximately 111km at the equator. Therefore if we
#' # travel 111km due east from the coordinate {0,0} we can verify that we have
#' # moved approximately 1 degree longitude and zero degrees latitude
#' bearing_to_lonlat(0, 0, 90, 111)

bearing_to_lonlat <- function(origin_lon, origin_lat, bearing, gc_dist) {
  
  # convert origin_lat, origin_lon and bearing from degrees to radians
  origin_lat <- origin_lat*2*pi/360
  origin_lon <- origin_lon*2*pi/360
  bearing <- bearing*2*pi/360
  
  # calculate new lat/lon using great circle distance
  earth_rad <- 6371
  new_lat <- asin(sin(origin_lat)*cos(gc_dist/earth_rad) + cos(origin_lat)*sin(gc_dist/earth_rad)*cos(bearing))
  new_lon <- origin_lon + atan2(sin(bearing)*sin(gc_dist/earth_rad)*cos(origin_lat), cos(gc_dist/earth_rad)-sin(origin_lat)*sin(new_lat))
  
  # convert new_lat and new_lon from radians to degrees
  new_lat <- new_lat*360/(2*pi)
  new_lon <- new_lon*360/(2*pi)
  
  return(list(longitude = new_lon,
              latitude = new_lat))
}

#------------------------------------------------
#' @title Calculate great circle distance and bearing between coordinates
#'
#' @description Calculate great circle distance and bearing between spatial
#'   coordinates.
#'
#' @param origin_lon The origin longitude
#' @param origin_lat The origin latitude
#' @param dest_lon The destination longitude
#' @param dest_lat The destination latitude
#'
#' @export
#' @examples
#' # one degree longitude should equal approximately 111km at the equator
#' lonlat_to_bearing(0, 0, 1, 0)

lonlat_to_bearing <- function(origin_lon, origin_lat, dest_lon, dest_lat) {
  
  # convert input arguments to radians
  origin_lon <- origin_lon*2*pi/360
  origin_lat <- origin_lat*2*pi/360
  dest_lon <- dest_lon*2*pi/360
  dest_lat <- dest_lat*2*pi/360
  
  delta_lon <- dest_lon - origin_lon
  
  # calculate bearing
  bearing <- atan2(sin(delta_lon)*cos(dest_lat), cos(origin_lat)*sin(dest_lat)-sin(origin_lat)*cos(dest_lat)*cos(delta_lon))
  
  # calculate great circle angle. Use temporary variable to avoid acos(>1) or 
  # acos(<0), which can happen due to underflow issues
  tmp <- sin(origin_lat)*sin(dest_lat) + cos(origin_lat)*cos(dest_lat)*cos(delta_lon)
  tmp <- ifelse(tmp > 1, 1, tmp)
  tmp <- ifelse(tmp < 0, 0, tmp)
  gc_angle <- acos(tmp)
  
  # convert bearing from radians to degrees measured clockwise from due north,
  # and convert gc_angle to great circle distance via radius of earth (km)
  bearing <- bearing*360/(2*pi)
  bearing <- (bearing+360)%%360
  earth_rad <- 6371
  gc_dist <- earth_rad*gc_angle
  
  # return list
  ret <-list(bearing = bearing,
             gc_dist = gc_dist)
  return(ret)
}

#------------------------------------------------
#' @title Convert lat/lon to cartesian coordinates
#'
#' @description Convert lat/lon coordinates to cartesian coordinates by first
#'   calculating great circle distance and bearing and then mapping these
#'   coordinates into cartesian space. This mapping is relative to the point
#'   {centre_lat, centre_lon}, which should be roughly at the midpoint of the
#'   observed data.
#'
#' @param centre_lon The centre longitude
#' @param centre_lat The centre latitude
#' @param data_lon The data longitude
#' @param data_lat The data latitude
#'
#' @export
#' @examples
#' # TODO

lonlat_to_cartesian <- function(centre_lon, centre_lat, data_lon, data_lat) {
  
  # calculate bearing and great circle distance of data relative to centre
  data_trans <- lonlat_to_bearing(centre_lon, centre_lat, data_lon, data_lat)
  
  # use bearing and distance to calculate cartesian coordinates
  theta <- data_trans$bearing*2*pi/360
  d <- data_trans$gc_dist
  data_x <- d*sin(theta)
  data_y <- d*cos(theta)
  
  # return list
  ret <- list(x = data_x,
              y = data_y)
  return(ret)
}

#------------------------------------------------
# Scaled Student's t distribution. Used in kernel density smoothing.
#' @noRd
dts <- function(x, df = 3, scale = 1, log = FALSE) {
  ret <- lgamma((df+1)/2) - lgamma(df/2) - 0.5*log(pi*df*scale^2) - ((df+1)/2)*log(1 + x^2/(df*scale^2))
  if (!log) {
    ret <- exp(ret)
  }
  return(ret)
}

#------------------------------------------------
# Bin values in two dimensions
#' @noRd
bin2D <- function(x, y, x_breaks, y_breaks) {
  
  # get number of breaks in each dimension
  nx <- length(x_breaks)
  ny <- length(y_breaks)
  
  # create table of binned values
  tab1 <- table(findInterval(x, x_breaks), findInterval(y, y_breaks))
  
  # convert to dataframe and force numeric
  df1 <- as.data.frame(tab1, stringsAsFactors = FALSE)
  names(df1) <- c("x", "y", "count")
  df1$x <- as.numeric(df1$x)
  df1$y <- as.numeric(df1$y)
  
  # subset to within breaks range
  df2 <- subset(df1, x > 0 & x < nx & y > 0 & y < ny)
  
  # fill in matrix
  mat1 <- matrix(0, ny-1, nx-1)
  mat1[cbind(df2$y, df2$x)] <- df2$count
  
  # calculate cell midpoints
  x_mids <- (x_breaks[-1] + x_breaks[-nx])/2
  y_mids <- (y_breaks[-1] + y_breaks[-ny])/2
  
  # return output as list
  ret <- list(x_mids = x_mids,
              y_mids = y_mids,
              z = mat1)
  return(ret)
}

#------------------------------------------------
#' Produce a smooth surface using 2D kernel density smoothing
#'
#' Takes lon/lat coordinates, bins in two dimensions and smooths using kernel 
#' density smoothing. Kernel densities are computed using the fast Fourier 
#' transform method, which is many times faster than simple summation when using
#' a large number of points. Each Kernel is student's-t distributed and scaled
#' by the bandwidth lambda. If lambda is set to \code{NULL} then the optimal
#' value of lambda is chosen automatically using the leave-one-out maximum
#' likelihood method.
#'
#' @param longitude longitude of input points
#' @param latitude latitude of input points
#' @param breaks_lon positions of longitude breaks
#' @param breaks_lat positions of latitude breaks
#' @param lambda bandwidth to use in posterior smoothing. If NULL then optimal 
#'   bandwidth is chosen automatically by maximum-likelihood
#' @param nu degrees of freedom of student's-t kernel
#'
#' @references Barnard, Etienne. "Maximum leave-one-out likelihood for kernel
#'   density estimation." Proceedings of the Twenty-First Annual Symposium of
#'   the Pattern Recognition Association of South Africa. 2010
#' @export
#' @examples
#' # TODO

kernel_smooth <- function(longitude, latitude, breaks_lon, breaks_lat, lambda = NULL, nu = 3) {
  
  # check inputs
  assert_numeric(longitude)
  assert_numeric(latitude)
  assert_same_length(longitude, latitude)
  assert_numeric(breaks_lon)
  assert_numeric(breaks_lat)
  if (!is.null(lambda)) {
    assert_single_pos(lambda, zero_allowed = FALSE)
  }
  assert_single_pos(nu, zero_allowed = FALSE)
  
  # get properties of cells in each dimension
  cells_lon <- length(breaks_lon) - 1
  cells_lat <- length(breaks_lat) - 1
  centre_lon <- mean(breaks_lon)
  centre_lat <- mean(breaks_lat)
  cellSize_lon <- diff(breaks_lon[1:2])
  cellSize_lat <- diff(breaks_lat[1:2])
  
  # bin lon/lat values in two dimensions and check that at least one value in
  # chosen region
  surface_raw <- bin2D(longitude, latitude, breaks_lon, breaks_lat)$z
  if (all(surface_raw == 0)) {
    stop('chosen lat/long window contains no posterior draws')
  }
  
  # temporarily add guard rail to surface to avoid Fourier series bleeding round
  # edges
  rail_size_lon <- cells_lon
  rail_size_lat <- cells_lat
  rail_mat_lon <- matrix(0, cells_lat, rail_size_lon)
  rail_mat_lat <- matrix(0, rail_size_lat, cells_lon + 2*rail_size_lon)
  
  surface_normalised <- surface_raw/sum(surface_raw)
  surface_normalised <- cbind(rail_mat_lon, surface_normalised, rail_mat_lon)
  surface_normalised <- rbind(rail_mat_lat, surface_normalised, rail_mat_lat)
  
  # calculate Fourier transform of posterior surface
  f1 = fftw2d(surface_normalised)
  
  # calculate x and y size of one cell in cartesian space. Because of
  # transformation, this size will technically be different for each cell, but
  # use centre of space to get a middling value
  cellSize_trans <- lonlat_to_cartesian(centre_lon, centre_lat, centre_lon + cellSize_lon, centre_lat + cellSize_lat)
  cellSize_trans_lon <- cellSize_trans$x
  cellSize_trans_lat <- cellSize_trans$y
  
  # produce surface over which kernel will be calculated. This surface wraps
  # around in both x and y (i.e. the kernel is actually defined over a torus)
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
  logLike <- -Inf
  for (i in 1:length(lambda_vec)) {
    
    # calculate Fourier transform of kernel
    lambda_this <- lambda_vec[i]
    kernel <- dts(kernel_s_mat, df=3, scale=lambda_this)
    f2 = fftw2d(kernel)
    
    # combine Fourier transformed surfaces and take inverse. f4 will ultimately
    # become the main surface of interest.
    f3 = f1*f2
    f4 = Re(fftw2d(f3, inverse = T))/length(surface_normalised)
    
    # subtract from f4 the probability density of each point measured from
    # itself. In other words, move towards a leave-one-out kernel density method
    f5 <- f4 - surface_normalised*dts(0, df = nu, scale = lambda_this)
    f5[f5<0] <- 0
    f5 <- f5/sum(f4)
    
    # calculate leave-one-out log-likelihood at each point on surface
    f6 <- surface_normalised*log(f5)
    
    # break if total log-likelihood is at a local maximum
    if (sum(f6, na.rm = T) < logLike) {
      break()
    }
    
    # otherwise update logLike
    logLike <- sum(f6,na.rm=T)
  }
  
  # report chosen value of lambda
  #if (is.null(lambda)) {
  #  message(sprintf("maximum likelihood lambda = %s", round(lambda_this,3)))
  #}
  
  # remove guard rail
  f4 <- f4[,(rail_size_lon+1):(ncol(f4)-rail_size_lon)]
  f4 <- f4[(rail_size_lat+1):(nrow(f4)-rail_size_lat),]
  
  # return surface
  return(f4)
}

#------------------------------------------------
#' @title Get ESS
#'
#' @description Returns effective sample size (ESS) of chosen model run.
#'
#' @param project an RgeoProfile project, as produced by the function 
#'   \code{rgeoprofile_project()}
#' @param K get ESS for this value of K
#'
#' @export
#' @examples
#' # TODO

get_ESS <- function(project, K = NULL) {
  
  # check inputs
  assert_custom_class(project, "rgeoprofile_project")
  if (!is.null(K)) {
    assert_single_pos_int(K, zero_allowed = FALSE)
  }
  
  # get active set and check non-zero
  s <- project$active_set
  if (s == 0) {
    stop("no active parameter set")
  }
  
  # set default K to first value with output
  null_output <- mapply(function(x) {is.null(x$summary$ESS)}, project$output$single_set[[s]]$single_K)
  if (all(null_output)) {
    stop("no ESS output for active parameter set")
  }
  if (is.null(K)) {
    K <- which(!null_output)[1]
    message(sprintf("using K = %s by default", K))
  }
  
  # check output exists for chosen K
  ESS <- project$output$single_set[[s]]$single_K[[K]]$summary$ESS
  if (is.null(ESS)) {
    stop(sprintf("no ESS output for K = %s of active set", K))
  }
  
  return(ESS)
}


##########################################################################################################
# MISC CLASSES

#------------------------------------------------
#' @title TODO
#'
#' @description custom print function for rgeoprofile_simdata.
#'
#' @param x TODO
#' @param ... TODO
#'
#' @export

print.rgeoprofile_simdata <- function(x, ...) {
  print(unclass(x))
  invisible(x)
}

#------------------------------------------------
#' @title TODO
#'
#' @description custom print function for rgeoprofile_qmatrix.
#'
#' @param x TODO
#' @param ... TODO
#'
#' @export

print.rgeoprofile_qmatrix <- function(x, ...) {
  print(unclass(x))
  invisible(x)
}
