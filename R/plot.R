
#------------------------------------------------
# red-to-blue colours
#' @noRd
col_hotcold <- function(n = 6) {
  raw_cols <- c("#D73027", "#FC8D59", "#FEE090", "#E0F3F8", "#91BFDB", "#4575B4")
  my_pal <- colorRampPalette(raw_cols)
  return(my_pal(n))
}

#------------------------------------------------
# blue-to-red colours. Full credit to tim.colors from the fields package, from 
# which these colours derive. Copied rather than including the fields package to
# avoid dependency on another package for the sake of a single colour scheme.
#' @noRd
col_tim <- function(n = 10) {
  raw_cols <- c("#00008F", "#00009F", "#0000AF", "#0000BF", 
                "#0000CF", "#0000DF", "#0000EF", "#0000FF", "#0010FF", 
                "#0020FF", "#0030FF", "#0040FF", "#0050FF", "#0060FF", 
                "#0070FF", "#0080FF", "#008FFF", "#009FFF", "#00AFFF", 
                "#00BFFF", "#00CFFF", "#00DFFF", "#00EFFF", "#00FFFF", 
                "#10FFEF", "#20FFDF", "#30FFCF", "#40FFBF", "#50FFAF", 
                "#60FF9F", "#70FF8F", "#80FF80", "#8FFF70", "#9FFF60", 
                "#AFFF50", "#BFFF40", "#CFFF30", "#DFFF20", "#EFFF10", 
                "#FFFF00", "#FFEF00", "#FFDF00", "#FFCF00", "#FFBF00", 
                "#FFAF00", "#FF9F00", "#FF8F00", "#FF8000", "#FF7000", 
                "#FF6000", "#FF5000", "#FF4000", "#FF3000", "#FF2000", 
                "#FF1000", "#FF0000", "#EF0000", "#DF0000", "#CF0000", 
                "#BF0000", "#AF0000", "#9F0000", "#8F0000", "#800000")
  my_pal <- colorRampPalette(raw_cols)
  return(my_pal(n))
}

#------------------------------------------------
#' @title Expand series of colours by interpolation
#'
#' @description Expand a series of colours by interpolation to produce any 
#'   number of colours from a given series. The pattern of interpolation is 
#'   designed so that (n+1)th value contains the nth value plus one more colour,
#'   rather than being a completely different series. For example, running
#'   \code{more_colours(5)} and \code{more_colours(4)}, the first 4 colours will
#'   be shared between the two series.
#'
#' @param n how many colours to return
#' @param raw_cols vector of colours to interpolate
#'
#' @export

more_colours <- function(n = 5, raw_cols = col_hotcold()) {
  
  # check inputs
  assert_single_pos_int(n, zero_allowed = FALSE)
  assert_string(raw_cols)
  assert_vector(raw_cols)
  
  # generate colour palette from raw colours
  my_palette <- colorRampPalette(raw_cols)
  
  # simple case if n small
  if (n <= 2) {
    return(my_palette(3)[1:n])
  }
  
  # interpolate colours by repeatedly splitting the [0,1] interval until we have
  # enough values. n_steps is the number of times we have to do this. n_breaks
  # is the number of breaks for each step
  n_steps <- ceiling(log(n-1)/log(2))
  n_breaks <- 2^(1:n_steps) + 1
  
  # split the [0,1] interval this many times and drop duplicated values
  s <- unlist(mapply(function(x) seq(0,1,l=x), n_breaks, SIMPLIFY = FALSE))
  s <- s[!duplicated(s)]
  
  # convert s to integer index
  w <- match(s, seq(0,1,l = n_breaks[n_steps]))
  w <- w[1:n]
  
  # get final colours
  all_cols <- my_palette(n_breaks[n_steps])
  ret <- all_cols[w]
  
  return(ret)
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
  assert_in(axis_type, 1:2)
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
  
  # get properties
  rungs <- nrow(loglike_intervals)
  
  # produce plot with different axis options
  plot1 <- ggplot(loglike_intervals) + theme_bw()
  if (axis_type == 1) {
    x_vec <- as.factor(1:rungs)
    plot1 <- plot1 + geom_segment(aes_(x = ~x_vec, y = ~Q2.5, xend = ~x_vec, yend = ~Q97.5))
    plot1 <- plot1 + geom_point(aes_(x = ~x_vec, y = ~Q50))
    plot1 <- plot1 + xlab("rung") + ylab("log-likelihood")
    
  } else if (axis_type == 2) {
    x_vec <- (1:rungs)/rungs
    plot1 <- plot1 + geom_segment(aes_(x = ~x_vec, y = ~Q2.5, xend = ~x_vec, yend = ~Q97.5))
    plot1 <- plot1 + geom_point(aes_(x = ~x_vec, y = ~Q50))
    plot1 <- plot1 + xlab(parse(text = "beta")) + ylab("log-likelihood")
    plot1 <- plot1 + coord_cartesian(xlim = c(0,1))
    
  } else {
    # TODO - make this option available if and when we implement temperature rungs
    #GTI_pow <- project$output$single_set[[s]]$single_K[[K]]$function_call$args$GTI_pow
    #x_vec <- ((1:rungs)/rungs)^GTI_pow
    #plot1 <- plot1 + geom_segment(aes_(x = ~x_vec, y = ~Q2.5, xend = ~x_vec, yend = ~Q97.5))
    #plot1 <- plot1 + geom_point(aes_(x = ~x_vec, y = ~Q50))
    #plot1 <- plot1 + xlab(parse(text = "beta^gamma")) + ylab("log-likelihood")
    #plot1 <- plot1 + coord_cartesian(xlim = c(0,1))
  }
  
  # optionally add central line
  if (connect_points) {
    plot1 <- plot1 + geom_line(aes(x = x_vec, y = loglike_intervals$Q50))
  }
  
  # optionally connect whiskers
  if (connect_whiskers) {
    plot1 <- plot1 + geom_line(aes(x = x_vec, y = loglike_intervals$Q2.5), linetype = "dotted") + geom_line(aes(x = x_vec, y = loglike_intervals$Q97.5), linetype = "dotted")
  }
  
  # return plot object
  return(plot1)
}

#------------------------------------------------
#' @title TODO
#'
#' @description default plot for class rgeoprofile_qmatrix.
#'
#' @param x TODO
#' @param y TODO
#' @param ... TODO
#'
#' @export

plot.rgeoprofile_qmatrix <- function(x, y, ...) {
  
  # get data into ggplot format
  m <- unclass(x)
  m <- m[!is.na(m[,1]), , drop = FALSE]
  n <- nrow(m)
  K <- ncol(m)
  df <- data.frame(sentinel_site = rep(1:n,each=K), k = as.factor(rep(1:K,times=n)), prob = as.vector(t(m)))
  
  # produce basic plot
  plot1 <- ggplot(df) + theme_empty()
  plot1 <- plot1 + geom_bar(aes_(x = ~sentinel_site, y = ~prob, fill = ~k), width = 1, stat = "identity")
  plot1 <- plot1 + scale_x_continuous(expand = c(0,0)) + scale_y_continuous(expand = c(0,0))
  plot1 <- plot1 + xlab("positive sentinel site") + ylab("probability")
  
  # add legends
  plot1 <- plot1 + scale_fill_manual(values = more_colours(K), name = "group")
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
  plot1 <- plot1 + xlab("positive sentinel site")
  
  # arrange in rows
  if (length(K) == 1) {
    plot1 <- plot1 + facet_wrap(~K, ncol = 1)
    plot1 <- plot1 + theme(strip.background = element_blank(), strip.text = element_blank())
    plot1 <- plot1 + ylab("probability")
  } else {
    plot1 <- plot1 + facet_wrap(~K, ncol = 1, strip.position = "left")
    plot1 <- plot1 + theme(strip.background = element_blank())
    plot1 <- plot1 + ylab("K")
  }
  
  # add legends
  plot1 <- plot1 + scale_fill_manual(values = more_colours(max(K)), name = "group")
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
  
  # get properties
  x_vec <- rownames(sigma_intervals)
  
  # produce plot
  plot1 <- ggplot(sigma_intervals) + theme_bw()
  plot1 <- plot1 + geom_segment(aes_(x = ~x_vec, y = ~Q2.5, xend = ~x_vec, yend = ~Q97.5))
  plot1 <- plot1 + geom_point(aes_(x = ~x_vec, y = ~Q50))
  plot1 <- plot1 + scale_y_continuous(limits = c(0, max(sigma_intervals$Q97.5)*1.1), expand = c(0,0))
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
  
  # produce plot
  plot1 <- ggplot(expected_popsize_intervals) + theme_bw()
  plot1 <- plot1 + geom_segment(aes_(x = "", y = ~Q2.5, xend = "", yend = ~Q97.5))
  plot1 <- plot1 + geom_point(aes_(x = 1, y = ~Q50))
  plot1 <- plot1 + scale_y_continuous(limits = c(0, max(expected_popsize_intervals$Q97.5)*1.1), expand = c(0,0))
  plot1 <- plot1 + xlab("") + ylab("expected population size")
  
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

#------------------------------------------------
#' @title Create dynamic map
#'
#' @description Create dynamic map
#'
#' @param map_type an index from 1 to 137 indicating the type of base map. The
#'   map types are taken from \code{leaflet::providers}. Defaults to "CartoDB"
#'
#' @export

plot_map <- function(map_type = 97) {
  
  # check inputs
  assert_in(map_type, 1:137, message = "map_type must be in 1:137")
  
  # produce plot
  myplot <- leaflet()
  myplot <-  addProviderTiles(myplot, leaflet::providers[[map_type]])
  
  # return plot object
  return(myplot)
}

#------------------------------------------------
#' @title Add sentinel sites to dynamic map
#'
#' @description Add sentinel sites to dynamic map
#'
#' @param myplot dynamic map produced by \code{plot_map()} function
#' @param project an RgeoProfile project, as produced by the function 
#'   \code{rgeoprofile_project()}
#' @param sentinel_radius the radius of sentinel sites. Taken from the active
#'   parameter set if unspecified
#' @param fill whether to fill circles
#' @param fill_colour colour of circle fill
#' @param fill_opacity fill opacity
#' @param border whether to add border to circles
#' @param border_colour colour of circle borders
#' @param border_weight thickness of circle borders
#' @param border_opacity opacity of circle borders
#'
#' @export

overlay_sentinels <- function(myplot,
                              project,
                              sentinel_radius = NULL,
                              fill = TRUE,
                              fill_colour = c(grey(0.5), "red"),
                              fill_opacity = 0.5,
                              border = FALSE,
                              border_colour = "black",
                              border_weight = 1,
                              border_opacity = 1.0) {
  
  # check inputs
  assert_custom_class(myplot, "leaflet")
  assert_custom_class(project, "rgeoprofile_project")
  if (!is.null(sentinel_radius)) {
    assert_single_pos(sentinel_radius)
  }
  assert_logical(fill)
  assert_vector(fill)
  assert_in(length(fill), c(1,2))
  if (length(fill) == 1) {
    fill <- rep(fill, 2)
  }
  assert_string(fill_colour)
  assert_vector(fill_colour)
  assert_in(length(fill_colour), c(1,2))
  if (length(fill_colour) == 1) {
    fill_colour <- rep(fill_colour, 2)
  }
  assert_single_pos(fill_opacity)
  assert_bounded(fill_opacity, 0, 1, inclusive_left = TRUE, inclusive_right = TRUE)
  assert_logical(border)
  assert_vector(border)
  assert_in(length(border), c(1,2))
  if (length(border) == 1) {
    border <- rep(border, 2)
  }
  assert_string(border_colour)
  assert_vector(border_colour)
  assert_in(length(border_colour), c(1,2))
  if (length(border_colour) == 1) {
    border_colour <- rep(border_colour, 2)
  }
  assert_single_pos(border_opacity)
  assert_bounded(border_opacity, 0, 1, inclusive_left = TRUE, inclusive_right = TRUE)
  
  # check for data
  df <- project$data
  if (is.null(df)) {
    stop("no data loaded")
  }
  
  # get sentinel radius from active parameter set by default
  if (is.null(sentinel_radius)) {
    message("getting sentinel radius from active parameter set:")
    
    # get active set and check non-zero
    s <- project$active_set
    if (s == 0) {
      stop("  no active parameter set")
    }
    
    # get sentinel radius
    sentinel_radius <- project$parameter_sets[[s]]$sentinel_radius
    message(sprintf("  sentinal radius = %skm", sentinel_radius))
  }
  
  # make circle attributes depend on counts
  n <- nrow(df)
  fill_vec <- rep(fill[1], n)
  fill_vec[df$counts > 0] <- fill[2]
  fill_colour_vec <- rep(fill_colour[1], n)
  fill_colour_vec[df$counts > 0] <- fill_colour[2]
  border_vec <- rep(border[1], n)
  border_vec[df$counts > 0] <- border[2]
  border_colour_vec <- rep(border_colour[1], n)
  border_colour_vec[df$counts > 0] <- border_colour[2]
  
  # overlay circles
  myplot <- addCircles(myplot, lng = df$longitude, lat = df$latitude,
                      radius = sentinel_radius*1e3,
                      fill = fill_vec, fillColor = fill_colour_vec, fillOpacity = fill_opacity,
                      stroke = border_vec, color = border_colour_vec,
                      opacity = border_opacity, weight = border_weight)
  
  # return plot object
  return(myplot)
}

#------------------------------------------------
#' @title Add points to dynamic map
#'
#' @description Add points to dynamic map
#'
#' @param myplot dynamic map produced by \code{plot_map()} function
#' @param lon longitude of points
#' @param lat latitude of points
#' @param col colour of points
#' @param size size of points
#' @param opacity opacity of points
#'
#' @export

overlay_points <- function(myplot, lon, lat, col = "black", size = 1, opacity = 1.0) {
  
  # check inputs
  assert_custom_class(myplot, "leaflet")
  assert_numeric(lon)
  assert_vector(lon)
  assert_numeric(lat)
  assert_vector(lat)
  assert_same_length(lon, lat)
  assert_single_string(col)
  assert_single_pos(size, zero_allowed = FALSE)
  
  # add circle markers
  myplot <- addCircleMarkers(myplot, lng = lon, lat = lat, radius = 2,
                             fillColor = col, stroke = FALSE, fillOpacity = opacity)
  
  # return plot object
  return(myplot)
}

#------------------------------------------------
#' @title Add geoprofile to dynamic map
#'
#' @description Add geoprofile to dynamic map
#'
#' @param myplot dynamic map produced by \code{plot_map()} function
#' @param project an RgeoProfile project, as produced by the function 
#'   \code{rgeoprofile_project()}
#' @param K which value of K to plot
#' @param source which source to plot. If NULL then plot combined surface
#' @param threshold what proportion of geoprofile to plot
#' @param col set of plotting colours
#' @param opacity opacity of geoprofile (that is not invisible due to being
#'   below threshold)
#'
#' @export

overlay_geoprofile <- function(myplot,
                               project,
                               K = NULL,
                               source = NULL,
                               threshold = 0.1,
                               col = col_hotcold(),
                               opacity = 0.8) {
  
  # check inputs
  assert_custom_class(myplot, "leaflet")
  assert_custom_class(project, "rgeoprofile_project")
  if (!is.null(source)) {
    assert_single_pos_int(source, zero_allowed = FALSE)
  }

  # get active set and check non-zero
  s <- project$active_set
  if (s == 0) {
    stop("  no active parameter set")
  }
  
  # set default K to first value with output
  null_output <- mapply(function(x) {is.null(x$summary$geoprofile)}, project$output$single_set[[s]]$single_K)
  if (all(null_output)) {
    stop("no geoprofile output for active parameter set")
  }
  if (is.null(K)) {
    K <- which(!null_output)[1]
    message(sprintf("using K = %s by default", K))
  }
  
  # check output exists for chosen K
  geoprofile <- project$output$single_set[[s]]$single_K[[K]]$summary$geoprofile
  if (is.null(geoprofile)) {
    stop(sprintf("no geoprofile output for K = %s of active set", K))
  }
  
  # choose which surface to plot
  if (is.null(source)) {
    source_plot <- "combined"
  } else {
    assert_leq(source, K)
    source_plot <- paste0("source", source)
  }
  
  # extract geoprofile into matrix
  gp <- t(matrix(geoprofile[[source_plot]], nrow = length(unique(geoprofile$lon))))
  gp[gp > threshold*100] <- NA
  
  # convert to raster
  r <- flip(raster(gp), direction = 2)
  
  # set extents and projection
  min_lon <- project$parameter_sets[[s]]$min_lon
  max_lon <- project$parameter_sets[[s]]$max_lon
  min_lat <- project$parameter_sets[[s]]$min_lat
  max_lat <- project$parameter_sets[[s]]$max_lat
  r <- setExtent(r, extent(min_lon, max_lon, min_lat, max_lat))
  crs(r) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" 
  
  # overlay raster
  myplot <- addRasterImage(myplot, x = r, colors = col, opacity = opacity)
  
  # add bounding rect
  myplot <- addRectangles(myplot, min_lon, min_lat, max_lon, max_lat,
                          fill = FALSE, weight = 2, color = grey(0.2))
  
  # return plot object
  return(myplot)
}

#------------------------------------------------
#' @title Add posterior probability surface to dynamic map
#'
#' @description Add posterior probability surface to dynamic map
#'
#' @param myplot dynamic map produced by \code{plot_map()} function
#' @param project an RgeoProfile project, as produced by the function 
#'   \code{rgeoprofile_project()}
#' @param K which value of K to plot
#' @param source which source to plot. If NULL then plot combined surface
#' @param threshold what proportion of posterior probability surface to plot
#' @param col set of plotting colours
#' @param opacity opacity of posterior probability surface (that is not
#'   invisible due to being below threshold)
#'
#' @export

overlay_surface <- function(myplot,
                            project,
                            K = NULL,
                            source = NULL,
                            threshold = 0.1,
                            col = col_hotcold(),
                            opacity = 0.8) {
  
  # check inputs
  assert_custom_class(myplot, "leaflet")
  assert_custom_class(project, "rgeoprofile_project")
  if (!is.null(source)) {
    assert_single_pos_int(source, zero_allowed = FALSE)
  }
  
  # get active set and check non-zero
  s <- project$active_set
  if (s == 0) {
    stop("  no active parameter set")
  }
  
  # set default K to first value with output
  null_output <- mapply(function(x) {is.null(x$summary$prob_surface)}, project$output$single_set[[s]]$single_K)
  if (all(null_output)) {
    stop("no prob_surface output for active parameter set")
  }
  if (is.null(K)) {
    K <- which(!null_output)[1]
    message(sprintf("using K = %s by default", K))
  }
  
  # check output exists for chosen K
  prob_surface <- project$output$single_set[[s]]$single_K[[K]]$summary$prob_surface
  if (is.null(prob_surface)) {
    stop(sprintf("no prob_surface output for K = %s of active set", K))
  }
  
  # choose which surface to plot
  if (is.null(source)) {
    source_plot <- "combined"
  } else {
    assert_leq(source, K)
    source_plot <- paste0("source", source)
  }
  
  # extract surface into matrix
  x <- t(matrix(prob_surface[[source_plot]], nrow = length(unique(prob_surface$lon))))
  x_min <- min(x, na.rm = TRUE)
  x_max <- max(x, na.rm = TRUE)
  threshold_final <- x_min + (x_max-x_min)*threshold
  x[x < threshold_final] <- NA
  
  # convert to raster
  r <- flip(raster(x), direction = 2)
  
  # set extents and projection
  min_lon <- project$parameter_sets[[s]]$min_lon
  max_lon <- project$parameter_sets[[s]]$max_lon
  min_lat <- project$parameter_sets[[s]]$min_lat
  max_lat <- project$parameter_sets[[s]]$max_lat
  r <- setExtent(r, extent(min_lon, max_lon, min_lat, max_lat))
  crs(r) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" 
  
  # overlay raster
  myplot <- addRasterImage(myplot, x = r, colors = rev(col), opacity = opacity)
  
  # add bounding rect
  myplot <- addRectangles(myplot, min_lon, min_lat, max_lon, max_lat,
                          fill = FALSE, weight = 2, color = grey(0.2))
  
  # return plot object
  return(myplot)
}

#------------------------------------------------
#' @title Add sources to dynamic map
#'
#' @description Add sources to dynamic map
#'
#' @param myplot dynamic map produced by \code{plot_map()} function
#' @param lon longitude of sources
#' @param lat latitude of sources
#' @param icon_url what image to use for the icon
#' @param icon_width the width of the icon
#' @param icon_height the height of the icon
#' @param icon_anchor_x the coordinates of the "tip" of the icon (relative to
#'   its top left corner, i.e. the top left corner means \code{icon_anchor_x =
#'   0} and \code{icon_anchor_y = 0}), and the icon will be aligned so that this
#'   point is at the marker's geographical location
#' @param icon_anchor_y the coordinates of the "tip" of the icon (relative to
#'   its top left corner, i.e. the top left corner means \code{icon_anchor_x =
#'   0} and \code{icon_anchor_y = 0}), and the icon will be aligned so that this
#'   point is at the marker's geographical location
#'
#' @export

overlay_sources <- function(myplot,
                            lon,
                            lat,
                            icon_url = "http://simpleicon.com/wp-content/uploads/cross.png",
                            icon_width = 20,
                            icon_height = 20,
                            icon_anchor_x = 10,
                            icon_anchor_y = 10) {
  
  # check inputs
  assert_custom_class(myplot, "leaflet")
  assert_numeric(lon)
  assert_vector(lon)
  assert_numeric(lat)
  assert_vector(lat)
  assert_same_length(lon, lat)
  assert_single_string(icon_url)
  assert_single_pos_int(icon_width)
  assert_single_pos_int(icon_height)
  assert_single_pos_int(icon_anchor_x)
  assert_single_pos_int(icon_anchor_y)
  
  # create custom icon
  source_icon <- makeIcon(iconUrl = icon_url, iconWidth = icon_width, iconHeight = icon_height,
                          iconAnchorX = icon_anchor_x, iconAnchorY = icon_anchor_y)
  
  # add custom markers
  myplot <- addMarkers(myplot, lng = lon, lat = lat, icon = source_icon)
  
  # return plot object
  return(myplot)
}

#------------------------------------------------
#' @title Add pie charts to dynamic map
#'
#' @description Add pie charts to dynamic map
#'
#' @param myplot dynamic map produced by \code{plot_map()} function
#' @param project an RgeoProfile project, as produced by the function 
#'   \code{rgeoprofile_project()}
#' @param K which value of K to plot
#' @param min_size lower limit on the size of pie charts
#' @param max_size upper limit on the size of pie charts
#' @param col segment colours
#'
#' @export

overlay_piecharts <- function(myplot,
                              project,
                              K = NULL,
                              min_size = 10,
                              max_size = 30,
                              col = NULL) {
  
  # check inputs
  assert_custom_class(myplot, "leaflet")
  assert_custom_class(project, "rgeoprofile_project")
  
  # get active set and check non-zero
  s <- project$active_set
  if (s == 0) {
    stop("  no active parameter set")
  }
  
  # set default K to first value with output
  null_output <- mapply(function(x) {is.null(x$summary$qmatrix)}, project$output$single_set[[s]]$single_K)
  if (all(null_output)) {
    stop("no qmatrix output for active parameter set")
  }
  if (is.null(K)) {
    K <- which(!null_output)[1]
    message(sprintf("using K = %s by default", K))
  }
  
  # check output exists for chosen K
  qmatrix <- project$output$single_set[[s]]$single_K[[K]]$summary$qmatrix
  if (is.null(qmatrix)) {
    stop(sprintf("no qmatrix output for K = %s of active set", K))
  }
  
  # set default colours from K
  col <- define_default(col, col_hotcold(K))
  
  # check correct number of colours
  assert_length(col, K)
  
  # get data into ggplot format
  w <- which(!is.na(qmatrix[,1]))
  lon <- project$data$longitude[w]
  lat <- project$data$latitude[w]
  counts <- project$data$counts[w]
  pie_size <- min_size + counts/max(counts)*(max_size - min_size)
  df <- round(qmatrix[w,], digits = 3)
  
  # overlay pie charts
  myplot <- addMinicharts(myplot, lon, lat,
                          type = "pie",
                          chartdata = df, 
                          colorPalette = col, 
                          width = pie_size,
                          transitionTime = 20)
  
  return(myplot)
}
