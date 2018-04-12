
#------------------------------------------------
#' Example London crimes
#'
#' Dummy crime sites simulated using rDPM() and centred on Queen Mary University of London.
#'
#' @docType data
#'
#' @format an object of class \code{list}.
#'
#' @usage data(LondonExample_crimes)
#'
#' @examples
#' # plot(LondonExample_crimes$latitude, LondonExample_crimes$longitude, xlab = "lon", ylab = "lat")
#' 
"LondonExample_crimes"

#------------------------------------------------
#' Example London sources
#'
#' Dummy sources associated with the crimes in LondonExample_crimes.
#'
#' @docType data
#'
#' @format an object of class \code{list}.
#'
#' @usage data(LondonExample_sources)
#' 
#' @examples
#' plot(LondonExample_crimes$latitude, LondonExample_crimes$longitude, xlab = "lon", ylab = "lat")
#' points(LondonExample_sources$latitude, LondonExample_sources$longitude, pch = 15, col = "blue")
#' 
"LondonExample_sources"

#------------------------------------------------
#' Cholera data
#'
#' Locations of cases of cholera from the 1854 cholera outbreak made famous by John Snow.
#'
#' @docType data
#' 
#' @format an object of class \code{list}.
#'
#' @usage data(Cholera)
#' 
#' @examples
#' plot(Cholera$latitude, Cholera$longitude, xlab = "lon", ylab = "lat")
#' 
"Cholera"

#------------------------------------------------

#' Cholera sources
#'
#' Locations of the neighbourhood water pumps near the site of the 1854 cholera outbreak. The Broad Street pump was the source of the outbreak.
#'
#' @docType data
#'
#' @format an object of class \code{list}.
#'
#' @usage data(WaterPumps)
#' 
#' @examples
#' plot(Cholera$latitude, Cholera$longitude, xlab = "lon", ylab = "lat")
#' points(WaterPumps$latitude, WaterPumps$longitude,  pch = 15, col = "blue")
#' 
"WaterPumps"

