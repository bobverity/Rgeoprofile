
#------------------------------------------------
# default RgeoProfile colours
#' @noRd
default_colours <- function(K) {
  
  # generate palette and colours
  raw_cols <- c("#D73027", "#FC8D59", "#FEE090", "#E0F3F8", "#91BFDB", "#4575B4")
  my_palette <- colorRampPalette(raw_cols)
  
  # simple case if small K
  if (K <= 2) {
    return(my_palette(K))
  }
  
  # some logic to choose a palette size and sequence of colours that is
  # consistent across different values of K
  ncol <- 3
  while(ncol<K) {
    ncol <- ncol+(ncol-1)
  }
  dist_mat <- matrix(1:ncol, ncol, ncol)
  dist_mat <- abs(t(dist_mat)-dist_mat)
  x <- rep(FALSE, ncol)
  
  col_index <- 1
  for (i in 2:K) {
    x[col_index] <- TRUE
    s <- apply(dist_mat[which(x),,drop=FALSE], 2, min)
    next_index <- which.max(s)
    col_index <- c(col_index, next_index)
  }
  col_index
  ret <- my_palette(ncol)[col_index]
  
  return(ret)
}

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
  data_all_observed <- subset(x$record$data_all, observed > 0)
  data_all_unobserved <- subset(x$record$data_all, observed == 0)
  
  # plot raw unobserved points
  plot1 <- ggplot() + theme_bw()
  plot1 <- plot1 + geom_point(aes_(x = ~longitude, y = ~latitude, col = "data_unobserved"),
                              size = 0.5, data = data_all_unobserved)
  
  # overlay true source locations
  plot1 <- plot1 + geom_point(aes_(x = ~longitude, y = ~latitude),
                              shape = 8, size = 2, col = "blue", data = x$record$true_source)
  
  # overlay circles around sentinel sites
  n_nodes <- 20
  for (i in 1:nrow(x$data)) {
    sentinel_lon <- x$data$longitude
    sentinel_lat <- x$data$latitude
    circle_lonlat <- as.data.frame(bearing_to_lonlat(sentinel_lon[i], sentinel_lat[i],
                                                     seq(0, 360, l=n_nodes), x$record$sentinel_radius))
    plot1 <- plot1 + geom_polygon(aes_(x = ~longitude, y = ~latitude),
                                  col = "#FF000099", fill = NA, data = circle_lonlat)
  }
  
  # overlay raw observed points
  plot1 <- plot1 + geom_point(aes_(x = ~longitude, y = ~latitude, col = "data_observed"),
                              size = 0.5, data = data_all_observed)
  
  # overlay count numbers around sentinel sites
  x$data$count_text <- mapply(function(x) {ifelse(x == 0, "", x)}, x$data$counts)
  plot1 <- plot1 + geom_text(aes_(x = ~longitude, y = ~latitude, label = ~count_text),
                             col = "red", data = x$data)
  
  # titles, legends, scales etc.
  plot1 <- plot1 + scale_color_manual(values = c("data_unobserved" = grey(0.7), "data_observed" = grey(0)))
  plot1 <- plot1 + xlab("longitude") + ylab("latitude")
  plot1 <- plot1 + guides(color = FALSE)
  
  return(plot1)
}

#------------------------------------------------
# ggplot theme with minimal objects
#' @noRd
theme_empty <- function() {
  theme(axis.line = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_blank())
}

#------------------------------------------------
# Default plot for class rgeoprofile_loglike_intervals
#' @noRd
plot.rgeoprofile_loglike_intervals <- function(x, y, ...) {
  
  # get data into ggplot format
  df <- as.data.frame(unclass(x))
  n <- nrow(df)
  x_vec <- y
  
  # produce plot
  plot1 <- ggplot(df) + theme_bw()
  plot1 <- plot1 + geom_segment(aes_(x = ~x_vec, y = ~Q2.5, xend = ~x_vec, yend = ~Q97.5))
  plot1 <- plot1 + geom_point(aes_(x = ~x_vec, y = ~Q50))
  plot1 <- plot1 + xlab("rung") + ylab("log-likelihood")
  
  # return plot object
  return(plot1)
}

#------------------------------------------------
#' @title Plot loglikelihood 95\% credible intervals
#'   
#' @description Plot loglikelihood 95\% credible intervals of current active set
#'   
#' @param project an RgeoProfile project, as produced by the function 
#'   \code{rgeoprofile_project()}
#' @param K which value of K to produce the plot for
#' @param axis_type how to format the x-axis. 1 = integer rungs, 2 = values of
#'   beta, 3 = values of beta raised to the GTI power
#' @param connect_points whether to connect points in the middle of intervals
#' @param connect_whiskers whether to connect points at the ends of the whiskers
#'
#' @export

plot_loglike <- function(project, K = NULL, axis_type = 1, connect_points = FALSE, connect_whiskers = FALSE) {
  
  # check inputs
  assert_custom_class(project, "rgeoprofile_project")
  if (!is.null(K)) {
    assert_single_pos_int(K, zero_allowed = FALSE)
  }
  assert_in(axis_type, 1:3)
  assert_single_logical(connect_points)
  assert_single_logical(connect_whiskers)
  
  # get active set and check non-zero
  s <- project$active_set
  if (s == 0) {
    stop("no active parameter set")
  }
  
  # set default K to first value with output
  null_output <- mapply(function(x) {is.null(x$summary$loglike_intervals)}, project$output$single_set[[s]]$single_K)
  if (all(null_output)) {
    stop("no loglike_intervals output for active parameter set")
  }
  if (is.null(K)) {
    K <- which(!null_output)[1]
    message(sprintf("using K = %s by default", K))
  }
  
  # check output exists for chosen K
  loglike_intervals <- project$output$single_set[[s]]$single_K[[K]]$summary$loglike_intervals
  if (is.null(loglike_intervals)) {
    stop(sprintf("no loglike_intervals output for K = %s of active set", K))
  }
  
  # produce plot with different axis options
  rungs <- nrow(loglike_intervals)
  if (axis_type == 1) {
    x_vec <- 1:rungs
    plot1 <- plot(loglike_intervals, as.factor(x_vec))
    
  } else if (axis_type == 2) {
    x_vec <- (1:rungs)/rungs
    plot1 <- plot(loglike_intervals, x_vec)
    plot1 <- plot1 + xlab(parse(text = "beta"))
    plot1 <- plot1 + coord_cartesian(xlim = c(0,1))
    
  } else {
    GTI_pow <- project$output$single_set[[s]]$single_K[[K]]$function_call$args$GTI_pow
    x_vec <- ((1:rungs)/rungs)^GTI_pow
    plot1 <- plot(loglike_intervals, x_vec)
    plot1 <- plot1 + xlab(parse(text = "beta^gamma"))
    plot1 <- plot1 + coord_cartesian(xlim = c(0,1))
  }
  
  # optionally add central line
  if (connect_points) {
    df <- as.data.frame(unclass(loglike_intervals))
    plot1 <- plot1 + geom_line(aes(x = x_vec, y = df$Q50))
  }
  
  # optionally connect whiskers
  if (connect_whiskers) {
    df <- as.data.frame(unclass(loglike_intervals))
    plot1 <- plot1 + geom_line(aes(x = x_vec, y = df$Q2.5), linetype = "dotted") + geom_line(aes(x = x_vec, y = df$Q97.5), linetype = "dotted")
  }
  
  # return plot object
  return(plot1)
}

#------------------------------------------------
#' @title Plot raw source locations
#'
#' @description Plot raw source locations.
#'
#' @details TODO
#'
#' @param project an RgeoProfile project, as produced by the function 
#'   \code{rgeoprofile_project()}
#' @param K which value of K to plot
#'
#' @export
#' @examples
#' # TODO

plot_source_raw <- function(project, K = NULL) {
  
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
  null_output <- mapply(function(x) {is.null(x$raw$source_lon)}, project$output$single_set[[s]]$single_K)
  if (all(null_output)) {
    stop("no source_lon output for active parameter set")
  }
  if (is.null(K)) {
    K <- which(!null_output)[1]
    message(sprintf("using K = %s by default", K))
  }
  
  # check output exists for chosen K
  source_lon <- project$output$single_set[[s]]$single_K[[K]]$raw$source_lon
  source_lat <- project$output$single_set[[s]]$single_K[[K]]$raw$source_lat
  if (is.null(source_lon) || is.null(source_lat)) {
    stop(sprintf("no source_lon/source_lat output for K = %s of active set", K))
  }
  
  # get raw sources into ggplot format
  df <- data.frame(source_lon = as.vector(source_lon),
                   source_lat = as.vector(source_lat),
                   source = rep(1:ncol(source_lon), each = nrow(source_lon)))
  
  # produce basic empty plot
  plot1 <- ggplot() + theme_bw()
  
  # overlay raw posterior source chains
  for (i in 1:ncol(source_lon)) {
    plot1 <- plot1 + geom_path(aes_(x = ~source_lon, y = ~source_lat, color = ~as.factor(source)), alpha = 0.2, data = df)
  }
  
  # overlay circles around sentinel sites
  n_nodes <- 20
  sentinel_radius <- project$parameter_sets[[s]]$sentinel_radius
  for (i in 1:nrow(project$data)) {
    sentinel_lon <- project$data$longitude
    sentinel_lat <- project$data$latitude
    circle_lonlat <- as.data.frame(bearing_to_lonlat(sentinel_lon[i], sentinel_lat[i],
                                                     seq(0, 360, l=n_nodes), sentinel_radius))
    plot1 <- plot1 + geom_polygon(aes_(x = ~longitude, y = ~latitude),
                                  col = grey(0.7), fill = NA, data = circle_lonlat)
  }
  
  # overlay count numbers around sentinel sites
  project$data$count_text <- mapply(function(x) {ifelse(x == 0, "", x)}, project$data$counts)
  plot1 <- plot1 + geom_text(aes_(x = ~longitude, y = ~latitude, label = ~count_text),
                             col = grey(0.7), data = project$data)
  
  # titles, legends, scales etc.
  plot1 <- plot1 + scale_color_manual(values = default_colours(K), name = "source")
  plot1 <- plot1 + xlab("longitude") + ylab("latitude")
  
  # return plot object
  return(plot1)
}

#------------------------------------------------
# Default plot for class rgeoprofile_qmatrix
#' @noRd
plot.rgeoprofile_qmatrix <- function(x, y, ...) {
  
  # get data into ggplot format
  m <- unclass(x)
  m <- m[!is.na(m[,1]),]
  n <- nrow(m)
  K <- ncol(m)
  df <- data.frame(ind = rep(1:n,each=K), k = as.factor(rep(1:K,times=n)), val = as.vector(t(m)))
  
  # produce basic plot
  plot1 <- ggplot(df) + theme_empty()
  plot1 <- plot1 + geom_bar(aes_(x = ~ind, y = ~val, fill = ~k), width = 1, stat = "identity")
  plot1 <- plot1 + scale_x_continuous(expand = c(0,0)) + scale_y_continuous(expand = c(0,0))
  plot1 <- plot1 + xlab("sample") + ylab("probability")
  
  # add legends
  plot1 <- plot1 + scale_fill_manual(values = default_colours(K), name = "group")
  plot1 <- plot1 + scale_colour_manual(values = "white")
  plot1 <- plot1 + guides(colour = FALSE)
  
  # add border
  plot1 <- plot1 + theme(panel.border = element_rect(colour = "black", size = 2, fill = NA))
  
  # return plot object
  return(plot1)
}

#------------------------------------------------
#' @title Posterior allocation plot
#'
#' @description Produce posterior allocation plot of current active set.
#'
#' @param project an RgeoProfile project, as produced by the function 
#'   \code{rgeoprofile_project()}
#' @param K which value of K to produce the plot for
#' @param divide_ind_on whether to add dividing lines between bars
#'
#' @export

plot_structure <- function(project, K = NULL, divide_ind_on = FALSE) {
  
  # check inputs
  assert_custom_class(project, "rgeoprofile_project")
  if (!is.null(K)) {
    assert_pos_int(K)
  }
  assert_single_logical(divide_ind_on)
  
  # get active set and check non-zero
  s <- project$active_set
  if (s == 0) {
    stop("no active parameter set")
  }
  
  # set default K to all values with output
  null_output <- mapply(function(x) {is.null(x$summary$qmatrix)}, project$output$single_set[[s]]$single_K)
  if (all(null_output)) {
    stop("no output for active parameter set")
  }
  K <- define_default(K, which(!null_output))
  
  # check output exists for chosen K
  qmatrix_list <- list()
  for (i in 1:length(K)) {
    qmatrix_list[[i]] <- project$output$single_set[[s]]$single_K[[K[i]]]$summary$qmatrix
    if (is.null(qmatrix_list[[i]])) {
      stop(sprintf("no qmatrix output for K = %s of active set", K[i]))
    }
  }
  
  # get data into ggplot format
  df <- NULL
  for (i in 1:length(K)) {
    m <- unclass(qmatrix_list[[i]])
    m <- m[!is.na(m[,1]), , drop = FALSE]
    n <- nrow(m)
    df <- rbind(df, data.frame(K = as.numeric(K[i]), ind = rep(1:n,each=K[i]), k = as.factor(rep(1:K[i],times=n)), val = as.vector(t(m))))
  }
  
  # produce basic plot
  plot1 <- ggplot(df) + theme_empty()
  plot1 <- plot1 + geom_bar(aes_(x = ~ind, y = ~val, fill = ~k), width = 1, stat = "identity")
  plot1 <- plot1 + scale_x_continuous(expand = c(0,0)) + scale_y_continuous(expand = c(0,0))
  
  # arrange in rows
  if (length(K) == 1) {
    plot1 <- plot1 + facet_wrap(~K, ncol = 1)
    plot1 <- plot1 + theme(strip.background = element_blank(), strip.text = element_blank())
    plot1 <- plot1 + xlab("sample") + ylab("probability")
  } else {
    plot1 <- plot1 + facet_wrap(~K, ncol = 1, strip.position = "left")
    plot1 <- plot1 + theme(strip.background = element_blank())
    plot1 <- plot1 + xlab("sample") + ylab("K")
  }
  
  # add legends
  plot1 <- plot1 + scale_fill_manual(values = default_colours(max(K)), name = "group")
  plot1 <- plot1 + scale_colour_manual(values = "white")
  plot1 <- plot1 + guides(colour = FALSE)
  
  # add border
  plot1 <- plot1 + theme(panel.border = element_rect(colour = "black", size = 2, fill = NA))
  
  # optionally add dividing lines
  if (divide_ind_on) {
    plot1 <- plot1 + geom_segment(aes_(x = ~x, y = ~y, xend = ~x, yend = ~y+1, col = "white"), size = 0.3, data = data.frame(x = 1:n-0.5, y = rep(0,n)))
  }
  
  return(plot1)
}

#------------------------------------------------
# Default plot for class rgeoprofile_sigma_intervals
#' @noRd
plot.rgeoprofile_sigma_intervals <- function(x, y, ...) {
  
  # get data into ggplot format
  df <- as.data.frame(unclass(x))
  n <- nrow(df)
  
  # produce plot
  plot1 <- ggplot(df) + theme_bw()
  plot1 <- plot1 + geom_segment(aes_(x = ~1:n, y = ~Q2.5, xend = ~1:n, yend = ~Q97.5))
  plot1 <- plot1 + geom_point(aes_(x = ~1:n, y = ~Q50))
  plot1 <- plot1 + scale_y_continuous(limits = c(0, max(df$Q97.5)*1.1), expand = c(0,0))
  plot1 <- plot1 + xlab("source") + ylab("sigma")
  
  # return plot object
  return(plot1)
}

#------------------------------------------------
#' @title Plot COI 95\% credible intervals
#'
#' @description Plot COI 95\% credible intervals of current active set
#'
#' @details TODO
#'
#' @param project an RgeoProfile project, as produced by the function 
#'   \code{rgeoprofile_project()}
#' @param K which value of K to plot
#'
#' @export
#' @examples
#' # TODO

plot_sigma <- function(project, K = NULL) {
  
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
  null_output <- mapply(function(x) {is.null(x$summary$sigma_intervals)}, project$output$single_set[[s]]$single_K)
  if (all(null_output)) {
    stop("no sigma_intervals output for active parameter set")
  }
  if (is.null(K)) {
    K <- which(!null_output)[1]
    message(sprintf("using K = %s by default", K))
  }
  
  # check output exists for chosen K
  sigma_intervals <- project$output$single_set[[s]]$single_K[[K]]$summary$sigma_intervals
  if (is.null(sigma_intervals)) {
    stop(sprintf("no sigma_intervals output for K = %s of active set", K))
  }
  
  # produce quantile plot
  plot1 <- plot(project$output$single_set[[s]]$single_K[[K]]$summary$sigma_intervals)
  
  # return plot object
  return(plot1)
}

#------------------------------------------------
# Default plot for class rgeoprofile_expected_popsize_intervals
#' @noRd
plot.rgeoprofile_expected_popsize_intervals <- function(x, y, ...) {
  
  # get data into ggplot format
  df <- data.frame(Q2.5 = x[1], Q50 = x[2], Q97.5 = x[3])
  
  # produce plot
  plot1 <- ggplot(df) + theme_bw()
  plot1 <- plot1 + geom_segment(aes_(x = "", y = ~Q2.5, xend = "", yend = ~Q97.5))
  plot1 <- plot1 + geom_point(aes_(x = 1, y = ~Q50))
  plot1 <- plot1 + scale_y_continuous(limits = c(0, max(df$Q97.5)*1.1), expand = c(0,0))
  plot1 <- plot1 + xlab("") + ylab("expected population size")
  
  # return plot object
  return(plot1)
}

#------------------------------------------------
#' @title Plot COI 95\% credible intervals
#'
#' @description Plot COI 95\% credible intervals of current active set
#'
#' @details TODO
#'
#' @param project an RgeoProfile project, as produced by the function 
#'   \code{rgeoprofile_project()}
#' @param K which value of K to plot
#'
#' @export
#' @examples
#' # TODO

plot_expected_popsize <- function(project, K = NULL) {
  
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
  null_output <- mapply(function(x) {is.null(x$summary$expected_popsize_intervals)}, project$output$single_set[[s]]$single_K)
  if (all(null_output)) {
    stop("no expected_popsize_intervals output for active parameter set")
  }
  if (is.null(K)) {
    K <- which(!null_output)[1]
    message(sprintf("using K = %s by default", K))
  }
  
  # check output exists for chosen K
  expected_popsize_intervals <- project$output$single_set[[s]]$single_K[[K]]$summary$expected_popsize_intervals
  if (is.null(expected_popsize_intervals)) {
    stop(sprintf("no expected_popsize_intervals output for K = %s of active set", K))
  }
  
  # produce quantile plot
  plot1 <- plot(project$output$single_set[[s]]$single_K[[K]]$summary$expected_popsize_intervals)
  
  # return plot object
  return(plot1)
}

#------------------------------------------------
#' @title Produce MCMC trace plot
#'   
#' @description Produce MCMC trace plot of the log-likelihood at each iteration.
#'   
#' @param project an RgeoProfile project, as produced by the function 
#'   \code{rgeoprofile_project()}
#' @param K which value of K to plot
#' @param rung which rung to plot. Defaults to the cold chain
#' @param col colour of the trace
#'   
#' @export

plot_trace <- function(project, K = NULL, rung = NULL, col = "black") {
  
  # check inputs
  assert_custom_class(project, "rgeoprofile_project")
  if (!is.null(K)) {
    assert_single_pos_int(K, zero_allowed = FALSE)
  }
  if (!is.null(rung)) {
    assert_single_pos_int(rung)
  }
  
  # get active set and check non-zero
  s <- project$active_set
  if (s == 0) {
    stop("no active parameter set")
  }
  
  # set default K to first value with output
  null_output <- mapply(function(x) {is.null(x$raw$loglike_sampling)}, project$output$single_set[[s]]$single_K)
  if (all(null_output)) {
    stop("no loglike_sampling output for active parameter set")
  }
  if (is.null(K)) {
    K <- which(!null_output)[1]
    message(sprintf("using K = %s by default", K))
  }
  
  # check output exists for chosen K
  loglike_sampling <- project$output$single_set[[s]]$single_K[[K]]$raw$loglike_sampling
  if (is.null(loglike_sampling)) {
    stop(sprintf("no loglike_sampling output for K = %s of active set", K))
  }
  
  # use cold rung by default
  rungs <- ncol(loglike_sampling)
  rung <- define_default(rung, rungs)
  assert_leq(rung, rungs)
  loglike <- as.vector(loglike_sampling[,rung])
  
  # get into ggplot format
  df <- data.frame(x = 1:length(loglike), y = loglike)
  
  # produce plot
  plot1 <- ggplot(df) + theme_bw() + ylab("log-likelihood")
  
  # complete plot
  plot1 <- plot1 + geom_line(aes_(x = ~x, y = ~y, colour = "col1"))
  plot1 <- plot1 + coord_cartesian(xlim = c(0,nrow(df)))
  plot1 <- plot1 + scale_x_continuous(expand = c(0,0))
  plot1 <- plot1 + scale_colour_manual(values = col)
  plot1 <- plot1 + guides(colour = FALSE)
  plot1 <- plot1 + xlab("iteration")
  
  # return plot object
  return(plot1)
}

#------------------------------------------------
#' @title Produce MCMC autocorrelation plot
#'
#' @description Produce MCMC autocorrelation plot of the log-likelihood
#'
#' @param project an RgeoProfile project, as produced by the function 
#'   \code{rgeoprofile_project()}
#' @param K which value of K to plot
#' @param rung which rung to plot. Defaults to the cold chain
#' @param col colour of the trace
#'
#' @export

plot_acf <- function(project, K = NULL, rung = NULL, col = "black") {
  
  # check inputs
  assert_custom_class(project, "rgeoprofile_project")
  if (!is.null(K)) {
    assert_single_pos_int(K, zero_allowed = FALSE)
  }
  if (!is.null(rung)) {
    assert_single_pos_int(rung)
  }
  
  # get active set and check non-zero
  s <- project$active_set
  if (s == 0) {
    stop("no active parameter set")
  }
  
  # set default K to first value with output
  null_output <- mapply(function(x) {is.null(x$raw$loglike_sampling)}, project$output$single_set[[s]]$single_K)
  if (all(null_output)) {
    stop("no loglike_sampling output for active parameter set")
  }
  if (is.null(K)) {
    K <- which(!null_output)[1]
    message(sprintf("using K = %s by default", K))
  }
  
  # check output exists for chosen K
  loglike_sampling <- project$output$single_set[[s]]$single_K[[K]]$raw$loglike_sampling
  if (is.null(loglike_sampling)) {
    stop(sprintf("no loglike_sampling output for K = %s of active set", K))
  }
  
  # use cold rung by default
  rungs <- ncol(loglike_sampling)
  rung <- define_default(rung, rungs)
  assert_leq(rung, rungs)
  loglike <- as.vector(loglike_sampling[,rung])
  
  # store variable to plot
  v <- loglike
  
  # get autocorrelation
  lag_max <- round(3*length(v)/effectiveSize(v))
  lag_max <- max(lag_max, 20)
  lag_max <- min(lag_max, length(v))
  
  # get into ggplot format
  a <- acf(v, lag.max = lag_max, plot = FALSE)
  acf <- as.vector(a$acf)
  df <- data.frame(lag = (1:length(acf))-1, ACF = acf)
  
  # produce plot
  plot1 <- ggplot(df) + theme_bw()
  plot1 <- plot1 + geom_segment(aes_(x = ~lag, y = 0, xend = ~lag, yend = ~ACF, colour = "col1"))
  plot1 <- plot1 + scale_colour_manual(values = col)
  plot1 <- plot1 + guides(colour = FALSE)
  plot1 <- plot1 + xlab("lag") + ylab("ACF")
  
  # return plot object
  return(plot1)
}

#------------------------------------------------
#' @title Produce MCMC density plot
#'
#' @description Produce MCMC density plot of the log-likelihood
#'
#' @param project an RgeoProfile project, as produced by the function 
#'   \code{rgeoprofile_project()}
#' @param K which value of K to plot
#' @param rung which rung to plot. Defaults to the cold chain
#' @param col colour of the trace
#'
#' @export

plot_density <- function(project, K = NULL, rung = NULL, col = "black") {
  
  # check inputs
  assert_custom_class(project, "rgeoprofile_project")
  if (!is.null(K)) {
    assert_single_pos_int(K, zero_allowed = FALSE)
  }
  if (!is.null(rung)) {
    assert_single_pos_int(rung)
  }
  
  # get active set and check non-zero
  s <- project$active_set
  if (s==0) {
    stop("no active parameter set")
  }
  
  # set default K to first value with output
  null_output <- mapply(function(x) {is.null(x$raw$loglike_sampling)}, project$output$single_set[[s]]$single_K)
  if (all(null_output)) {
    stop("no loglike_sampling output for active parameter set")
  }
  if (is.null(K)) {
    K <- which(!null_output)[1]
    message(sprintf("using K = %s by default", K))
  }
  
  # check output exists for chosen K
  loglike_sampling <- project$output$single_set[[s]]$single_K[[K]]$raw$loglike_sampling
  if (is.null(loglike_sampling)) {
    stop(sprintf("no loglike_sampling output for K = %s of active set", K))
  }
  
  # use cold rung by default
  rungs <- ncol(loglike_sampling)
  rung <- define_default(rung, rungs)
  assert_leq(rung, rungs)
  loglike <- as.vector(loglike_sampling[,rung])
  
  # get into ggplot format
  df <- data.frame(v = loglike)
  
  # produce plot
  plot1 <- ggplot(df) + theme_bw() + xlab("log-likelihood")
  
  # produce plot
  #plot1 <- ggplot(df) + theme_bw()
  plot1 <- plot1 + geom_histogram(aes_(x = ~v, y = ~..density.., fill = "col1"), bins = 50)
  plot1 <- plot1 + scale_fill_manual(values = col)
  plot1 <- plot1 + guides(fill = FALSE)
  plot1 <- plot1 + ylab("density")
  
  # return plot object
  return(plot1)
}

#------------------------------------------------
#' @title Produce diagnostic plots of log-likelihood
#'
#' @description Produce diagnostic plots of the log-likelihood.
#'
#' @param project an RgeoProfile project, as produced by the function 
#'   \code{rgeoprofile_project()}
#' @param K which value of K to plot
#' @param rung which rung to plot. Defaults to the cold chain
#' @param col colour of the trace
#'
#' @export

plot_loglike_dignostic <- function(project, K = NULL, rung = NULL, col = "black") {
  
  # check inputs
  assert_custom_class(project, "rgeoprofile_project")
  if (!is.null(K)) {
    assert_single_pos_int(K, zero_allowed = FALSE)
  }
  if (!is.null(rung)) {
    assert_single_pos_int(rung)
  }
  
  # get active set and check non-zero
  s <- project$active_set
  if (s == 0) {
    stop("no active parameter set")
  }
  
  # set default K to first value with output
  null_output <- mapply(function(x) {is.null(x$raw$loglike_sampling)}, project$output$single_set[[s]]$single_K)
  if (all(null_output)) {
    stop("no loglike_sampling output for active parameter set")
  }
  if (is.null(K)) {
    K <- which(!null_output)[1]
    message(sprintf("using K = %s by default", K))
  }
  
  # check output exists for chosen K
  loglike_sampling <- project$output$single_set[[s]]$single_K[[K]]$raw$loglike_sampling
  if (is.null(loglike_sampling)) {
    stop(sprintf("no loglike_sampling output for K = %s of active set", K))
  }
  
  # use cold rung by default
  rungs <- ncol(loglike_sampling)
  rung <- define_default(rung, rungs)
  assert_leq(rung, rungs)
  
  # produce individual diagnostic plots and add features
  plot1 <- plot_trace(project, K = K, rung = rung, col = col)
  plot1 <- plot1 + ggtitle("MCMC trace")
  
  plot2 <- plot_acf(project, K = K, rung = rung, col = col)
  plot2 <- plot2 + ggtitle("autocorrelation")
  
  plot3 <- plot_density(project, K = K, rung = rung, col = col)
  plot3 <- plot3 + ggtitle("density")
  
  # produce grid of plots
  ret <- grid.arrange(plot1, plot2, plot3, layout_matrix = rbind(c(1,1), c(2,3)))
}
