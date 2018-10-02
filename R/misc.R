
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
  if (i==max_i) {
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
  
  # check for exact equality of points
  if (origin_lon == dest_lon && origin_lat == dest_lat) {
    return(list(bearing = 0, gc_dist = 0))
  }
  
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
  
  return(list(bearing = bearing, gc_dist = gc_dist))
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
  
  # print raw list
  print(unclass(x))
  
  # return invisibly
  invisible(x)
}

#------------------------------------------------
# Overload print function for malecot_qmatrix
#' @noRd
print.rgeoprofile_qmatrix <- function(x, ...) {
  
  # print raw list
  print(unclass(x))
  
  # return invisibly
  invisible(x)
}
