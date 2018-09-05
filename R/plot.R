
#------------------------------------------------
#' @title TODO
#'
#' @description default plot for class rgeoprofile_simdata.
#'
#' @param x TODO
#' @param y TODO
#' @param ... TODO
#'
#' @export

plot.rgeoprofile_simdata <- function(x, y, ...) {
  
  # subset observed vs. unobserved data
  data_all_observed <- subset(x$data_all, observed > 0)
  data_all_unobserved <- subset(x$data_all, observed == 0)
  
  # plot raw unobserved points
  plot1 <- ggplot() + theme_bw()
  plot1 <- plot1 + geom_point(aes_(x = ~longitude, y = ~latitude, col = "data_unobserved"),
                              size = 0.5, data = data_all_unobserved)
  
  # overlay true source locations
  plot1 <- plot1 + geom_point(aes_(x = ~longitude, y = ~latitude),
                              shape = 8, size = 2, col = "blue", data = x$source)
  
  # overlay circles around sentinel sites
  n_nodes <- 20
  for (i in 1:nrow(x$data_observed)) {
    sentinel_lon <- x$data_observed$longitude
    sentinel_lat <- x$data_observed$latitude
    circle_lonlat <- as.data.frame(bearing_to_lonlat(sentinel_lon[i], sentinel_lat[i],
                                                     seq(0, 360, l=n_nodes), x$sentinel_radius))
    plot1 <- plot1 + geom_polygon(aes_(x = ~longitude, y = ~latitude),
                                  col = "#FF000099", fill = NA, data = circle_lonlat)
  }
  
  # overlay raw observed points
  plot1 <- plot1 + geom_point(aes_(x = ~longitude, y = ~latitude, col = "data_observed"),
                              size = 0.5, data = data_all_observed)
  
  # overlay count numbers around sentinel sites
  x$data_observed$count_text <- mapply(function(x) {ifelse(x == 0, "", x)}, x$data_observed$counts)
  plot1 <- plot1 + geom_text(aes_(x = ~longitude, y = ~latitude, label = ~count_text),
                             col = "red", data = x$data_observed)
  
  # titles, legends, scales etc.
  plot1 <- plot1 + scale_color_manual(values = c("data_unobserved" = grey(0.7), "data_observed" = grey(0)))
  plot1 <- plot1 + xlab("longitude") + ylab("latitude")
  plot1 <- plot1 + guides(color = FALSE)
  
  return(plot1)
}

