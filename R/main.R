
#------------------------------------------------
# The following commands ensure that package dependencies are listed in the
# NAMESPACE file.

# Rcpp          - allows C++ integration
# parallel      - running jobs in parallel
# coda          - MCMC type objects and methods
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
#' @import parallel
#' @import coda
#' @importFrom fftwtools fftw2d
#' @import ggplot2
#' @import gridExtra
#' @import RColorBrewer
#' @import leaflet
#' @import leaflet.minicharts
#' @import rgdal
#' @importFrom raster raster values setValues addLayer extract xmin xmax ymin ymax disaggregate flip crs<- setExtent extent extent<- rasterize projectRaster distance
#' @import viridis
#' @importFrom grDevices colorRampPalette grey
#' @import graphics
#' @import stats
#' @import utils
NULL

#------------------------------------------------
#' @title Check that RgeoProfile package has loaded successfully
#'
#' @description Simple function to check that RgeoProfile package has loaded 
#'   successfully. Prints "RgeoProfile loaded successfully!" if so.
#'
#' @export

check_rgeoprofile_loaded <- function() {
  message("RgeoProfile loaded successfully!")
}

#------------------------------------------------
#' @title Bind data to project
#'   
#' @description Load data into a \code{rgeoprofile_project} prior to analysis.
#'   Data must be formatted as a dataframe with samples in rows and loci in
#'   columns. If individuals are polyploid then multiple rows can be used per
#'   sample. Ploidy is allowed to vary between samples, and can be specified in
#'   multiple ways.
#'   
#' @param project an \code{rgeoprofile_project}, as produced by the function
#'   \code{rgeoprofile_project()}
#' @param df a dataframe with three columns named "longitude", "latitude",
#'   "counts" (in that order). Counts must be positive integers or zero
#' @param name optional name of the data set to aid in record keeping
#' @param check_delete_output whether to prompt the user before overwriting 
#'   existing data
#'
#' @export

bind_data <- function(project, df, name = NULL, check_delete_output = TRUE) {
  
  # check inputs
  assert_custom_class(project, "rgeoprofile_project")
  assert_dataframe(df)
  assert_ncol(df, 3)
  assert_eq(names(df), c("longitude", "latitude", "counts"))
  assert_pos_int(df$counts, zero_allowed = TRUE)
  if (!is.null(name)) {
    assert_single_string(name)
  }
  assert_single_logical(check_delete_output)
  
  # check before overwriting existing output
  if (project$active_set>0 && check_delete_output) {
    
    # ask before overwriting. On abort, return original project
    if (!user_yes_no("All existing output and parameter sets for this project will be lost. Continue? (Y/N): ")) {
      return(project)
    }
    
    # replace old project with fresh empty version
    project <- rgeoprofile_project()
  }
  
  # add data to project
  project$data <- df
  
  return(project)
}

#------------------------------------------------
#' @title Create new parameter set
#'   
#' @description Create a new parameter set within an \code{rgeoprofile_project}. The new 
#'   parameter set becomes the active set once created.
#'   
#' @param project an rgeoprofile_project, as produced by the function 
#'   \code{rgeoprofile_project()}
#' @param name an optional name for the parameter set
#' @param sentinel_radius the observation radius of sentinel sites
#' @param min_lon the minimum possible longitude of source locations
#' @param max_lon the maximum possible longitude of source locations
#' @param min_lat the minimum possible latitude of source locations
#' @param max_lat the maximum possible latitude of source locations
#' @param sigma_model set as \code{"single"} to assume the same dispersal
#'   distance for all sources, or \code{"independent"} to assume an
#'   independently drawn dispersal distance for each source
#' @param sigma_prior_mean the prior mean of the parameter sigma (km)
#' @param sigma_prior_sd the prior standard deviation of the parameter sigma 
#'   (km). Set to 0 to use a fixed value for sigma (fixed at
#'   \code{sigma_prior_mean})
#' @param expected_popsize_prior_mean the prior mean of the expected total
#'   population size
#' @param expected_popsize_prior_sd the prior standard deviation of the expected
#'   total population size. Set to 0 to use a fixed value (fixed at
#'   \code{expected_popsize_prior_mean}), or set to -1 to use an improper,
#'   infinitely diffuse prior
#' 
#' @export

new_set <- function(project,
                    name = "(no name)",
                    sentinel_radius = 0.2,
                    min_lon = -10,
                    max_lon = 10,
                    min_lat = -10,
                    max_lat = 10,
                    sigma_model = "single",
                    sigma_prior_mean = 1,
                    sigma_prior_sd = 1,
                    expected_popsize_prior_mean = 100,
                    expected_popsize_prior_sd = 10) {
  
  # check inputs
  assert_custom_class(project, "rgeoprofile_project")
  assert_single_pos(sentinel_radius, zero_allowed = FALSE)
  assert_single_numeric(min_lon)
  assert_single_numeric(max_lon)
  assert_gr(max_lon, min_lon)
  assert_single_numeric(min_lat)
  assert_single_numeric(max_lat)
  assert_gr(max_lat, min_lat)
  assert_in(sigma_model, c("single", "independent"))
  assert_single_pos(sigma_prior_mean, zero_allowed = FALSE)
  assert_single_pos(sigma_prior_sd, zero_allowed = TRUE)
  assert_single_pos(expected_popsize_prior_mean, zero_allowed = FALSE)
  if (expected_popsize_prior_sd != -1) {
    assert_single_pos(expected_popsize_prior_sd, zero_allowed = TRUE)
  }
  
  # count current parameter sets and add one
  s <- length(project$parameter_sets) + 1
  
  # make new set active
  project$active_set <- s
  
  # create new parameter set
  project$parameter_sets[[s]] <- list(name = name,
                                      sentinel_radius = sentinel_radius,
                                      min_lon = min_lon,
                                      max_lon = max_lon,
                                      min_lat = min_lat,
                                      max_lat = max_lat,
                                      sigma_model = sigma_model,
                                      sigma_prior_mean = sigma_prior_mean,
                                      sigma_prior_sd = sigma_prior_sd,
                                      expected_popsize_prior_mean = expected_popsize_prior_mean,
                                      expected_popsize_prior_sd = expected_popsize_prior_sd)
  
  # name parameter set
  names(project$parameter_sets)[s] <- paste0("set", s)
  
  # create new output at all_K level
  project$output$single_set[[s]] <- list(single_K = list(),
                                         all_K = list())
  
  # name new output
  names(project$output$single_set) <- paste0("set", 1:length(project$output$single_set))
  
  # return
  return(project)
}

#------------------------------------------------
#' @title Delete parameter set
#'   
#' @description Delete a given parameter set from an \code{rgeoprofile_project}.
#'   
#' @param project an rgeoprofile_project, as produced by the function 
#'   \code{rgeoprofile_project()}
#' @param set which set to delete. Defaults to the current active set
#' @param check_delete_output whether to prompt the user before deleting any
#'   existing output
#'   
#' @export

delete_set <- function(project, set = NULL, check_delete_output = TRUE) {
  
  # check inputs
  assert_custom_class(project, "rgeoprofile_project")
  assert_single_logical(check_delete_output)
  
  # set index to active_set by default
  set <- define_default(set, project$active_set)
  
  # further checks
  assert_single_pos_int(set, zero_allowed = FALSE)
  assert_leq(set, length(project$parameter_sets))
  
  # check before overwriting existing output
  if (project$active_set>0 & check_delete_output) {
    
    # ask before overwriting. On abort, return original project
    if (!user_yes_no(sprintf("Any existing output for set %s will be deleted. Continue? (Y/N): ", set))) {
      return(project)
    }
  }
  
  # drop chosen parameter set
  project$parameter_sets[[set]] <- NULL
  
  # drop chosen output
  project$output$single_set[[set]] <- NULL
  
  # make new final set active
  project$active_set <- length(project$parameter_sets)
  
  # return
  return(project)
}

#------------------------------------------------
#' @title Run main MCMC
#'   
#' @description Run the main RgeoProfile MCMC. Model parameters are taken from
#'   the current active parameter set, and MCMC parameters are passed in as 
#'   arguments. All output is stored within the project.
#'   
#' @details Both longitude and latitude values can be represented to a given 
#'   precision level using the arguments \code{precision_lon} and 
#'   \code{precision_lat} - for example, a precision of 0.01 means that values 
#'   are rounded to the nearest hundredth of a degree. This allows the use of 
#'   look-up tables for the likelihood calculation, which significantly speeds 
#'   up the MCMC. Set to 0 to use exact values (up to C++ "double" precision)
#'   rather than using look-up tables.
#'   
#' @param project an rgeoprofile_project, as produced by the function 
#'   \code{rgeoprofile_project()}
#' @param K the number of sources
#' @param precision_lon the level of precision at which longitudes are
#'   represented
#' @param precision_lat the level of precision at which latitudes are
#'   represented
#' @param burnin the number of burn-in iterations
#' @param samples the number of sampling iterations
#' @param auto_converge whether convergence should be assessed automatically 
#'   every \code{converge_test} iterations, leading to termination of the 
#'   burn-in phase. If \code{FALSE} then the full \code{burnin} iterations are 
#'   used
#' @param converge_test test for convergence every \code{convergence_test} 
#'   iterations if \code{auto_converge} is being used
#' @param solve_label_switching_on whether to implement the Stevens' solution to
#'   the label-switching problem. If turned off then Q-matrix output will no 
#'   longer be correct, although evidence estimates will be unaffected
#' @param cluster option to pass in a cluster environment (see package 
#'   "parallel")
#' @param pb_markdown whether to run progress bars in markdown mode, in which 
#'   case they are updated once at the end to avoid large amounts of output
#' @param store_raw whether to store raw MCMC output in addition to summary 
#'   output. Setting to FALSE can considerably reduce output size in memory
#' @param silent whether to suppress all console output
#' 
#' @export

run_mcmc <- function(project, K = 3, precision_lon = 1e-3, precision_lat = 1e-3, burnin = 1e2, samples = 1e3, auto_converge = TRUE, converge_test = ceiling(burnin/10), solve_label_switching_on = TRUE, cluster = NULL, pb_markdown = FALSE, store_raw = TRUE, silent = !is.null(cluster)) {
  
  # start timer
  t0 <- Sys.time()
  
  # check inputs
  assert_custom_class(project, "rgeoprofile_project")
  assert_pos_int(K, zero_allowed = FALSE)
  assert_single_pos(precision_lon, zero_allowed = TRUE)
  assert_single_pos(precision_lat, zero_allowed = TRUE)
  assert_single_pos_int(burnin, zero_allowed = FALSE)
  assert_single_pos_int(samples, zero_allowed = FALSE)
  assert_single_logical(auto_converge)
  assert_single_pos_int(converge_test, zero_allowed = FALSE)
  assert_single_logical(solve_label_switching_on)
  if (!is.null(cluster)) {
    assert_custom_class(project, "cluster")
  }
  assert_single_logical(pb_markdown)
  assert_single_logical(silent)
  
  # get active set
  s <- project$active_set
  if (s == 0) {
    stop("no active parameter set")
  }
  
  # get user-defined domain limits
  min_lon <- project$parameter_sets[[s]]$min_lon
  max_lon <- project$parameter_sets[[s]]$max_lon
  min_lat <- project$parameter_sets[[s]]$min_lat
  max_lat <- project$parameter_sets[[s]]$max_lat
  
  # get final domain limits used within model
  min_lon_model <- ceiling(min_lon/precision_lon)*precision_lon
  max_lon_model <- floor(max_lon/precision_lon)*precision_lon
  min_lat_model <- ceiling(min_lat/precision_lat)*precision_lat
  max_lat_model <- floor(max_lat/precision_lat)*precision_lat
  
  # ---------- create argument lists ----------
  
  # data list
  args_data <- list(longitude = project$data$longitude,
                    latitude = project$data$latitude,
                    counts = project$data$counts)
  
  # input arguments list
  args_inputs <- list(precision_lon = precision_lon,
                      precision_lat = precision_lat,
                      burnin = burnin,
                      samples = samples,
                      auto_converge = auto_converge,
                      converge_test = converge_test,
                      solve_label_switching_on = solve_label_switching_on,
                      pb_markdown = pb_markdown,
                      silent = silent)
  
  # combine model parameters list with input arguments
  args_model <- c(project$parameter_sets[[s]], args_inputs)
  args_model$sigma_model_numeric <- match(args_model$sigma_model, c("single", "independent"))
  args_model <- c(args_model, list(min_lon_model = min_lon_model,
                                   max_lon_model = max_lon_model,
                                   min_lat_model = min_lat_model,
                                   max_lat_model = max_lat_model))
  
  # R functions to pass to Rcpp
  args_functions <- list(test_convergence = test_convergence,
                         update_progress = update_progress)
  
  # define final argument list over all K
  parallel_args <- list()
  for (i in 1:length(K)) {
    
    # create progress bars
    pb_burnin <- txtProgressBar(min = 0, max = burnin, initial = NA, style = 3)
    pb_samples <- txtProgressBar(min = 0, max = samples, initial = NA, style = 3)
    args_progress <- list(pb_burnin = pb_burnin,
                          pb_samples = pb_samples)
    
    # incporporate arguments unique to this K
    args_model$K <- K[i]
    
    # create argument list
    parallel_args[[i]] <- list(args_data = args_data,
                               args_model = args_model,
                               args_functions = args_functions,
                               args_progress = args_progress)
  }
  
  # ---------- run MCMC ----------
  
  # split into parallel and serial implementations
  if (!is.null(cluster)) { # run in parallel
    clusterEvalQ(cluster, library(rmaverick))
    output_raw <- clusterApplyLB(cl = cluster, parallel_args, run_mcmc_cpp)
  } else { # run in serial
    output_raw <- lapply(parallel_args, run_mcmc_cpp)
  }
  
  #------------------------
  
  # begin processing results
  if (!silent) {
    cat("Processing results\n")
  }
  
  # loop through K
  ret <- list()
  all_converged <- TRUE
  for (i in 1:length(K)) {
    
    # create name lists
    rungs <- 1
    group_names <- paste0("group", 1:K[i])
    rung_names <- paste0("rung", 1:rungs)
    
    # ---------- raw mcmc results ----------
    
    # get loglikelihood in coda::mcmc format
    loglike_burnin <- mapply(function(x){mcmc(x)}, output_raw[[i]]$loglike_burnin)
    colnames(loglike_burnin) <- rung_names
    loglike_sampling <- mcmc(t(rcpp_to_mat(output_raw[[i]]$loglike_sampling)))
    colnames(loglike_sampling) <- rung_names
    
    # get source lon lat in coda::mcmc format
    full_source_lon <- mcmc(rcpp_to_mat(output_raw[[i]]$source_lon))
    colnames(full_source_lon) <- group_names
    full_source_lat <- mcmc(rcpp_to_mat(output_raw[[i]]$source_lat))
    colnames(full_source_lat) <- group_names
    
    # get sigma in coda::mcmc format
    full_sigma <- mcmc(rcpp_to_mat(output_raw[[i]]$sigma))
    if (args_model$sigma_model == "single") {
      full_sigma <- full_sigma[, 1, drop = FALSE]
      colnames(full_sigma) <- "all_groups"
    } else {
      colnames(full_sigma) <- group_names
    }
    
    # get expected_popsize in coda::mcmc format
    full_expected_popsize <- mcmc(output_raw[[i]]$expected_popsize)
    
    # ---------- summary results ----------
    
    # get 95% credible intervals over sampling loglikelihoods
    loglike_intervals <- as.data.frame(t(apply(loglike_sampling, 2, quantile_95)))
    
    # get 95% credible intervals over sigma
    sigma_intervals <- as.data.frame(t(apply(full_sigma, 2, quantile_95)))
    
    # get 95% credible intervals over expected_popsize
    expected_popsize_intervals <- as.data.frame(t(quantile_95(full_expected_popsize)))
    rownames(expected_popsize_intervals) <- "expected_popsize"
    
    # process Q-matrix
    qmatrix <- rcpp_to_mat(output_raw[[i]]$qmatrix)/samples
    qmatrix[project$data$counts == 0,] <- rep(NA, K[i])
    colnames(qmatrix) <- group_names
    class(qmatrix) <- "rgeoprofile_qmatrix"
    
    # create empty raster with correct properties
    cells_lon <- (max_lon_model - min_lon_model)/precision_lon
    cells_lat <- (max_lat_model - min_lat_model)/precision_lat
    raster_empty <- raster(xmn = min_lon_model,
                           xmx = max_lon_model,
                           ymn = min_lat_model,
                           ymx = max_lat_model,
                           nrow = cells_lat,
                           ncol = cells_lon)
    
    # get breaks for kernel smoothing
    breaks_lon <- seq(min_lon_model, max_lon_model, precision_lon)
    breaks_lat <- seq(min_lat_model, max_lat_model, precision_lat)
    
    # produce posterior probability surface rasters
    prob_surface_split <- raster()
    prob_surface_mat <- 0
    for (k in 1:K[i]) {
      
      # get prob_surface for this K by smoothing
      prob_surface_split_mat <- kernel_smooth(full_source_lon[,k],
                                              full_source_lat[,k],
                                              breaks_lon,
                                              breaks_lat)
      prob_surface_split_mat <- prob_surface_split_mat[nrow(prob_surface_split_mat):1,]
      prob_surface_split_mat <- prob_surface_split_mat/sum(prob_surface_split_mat)
      
      # add raster layer
      prob_surface_split_k <- setValues(raster_empty, prob_surface_split_mat)
      prob_surface_split <- addLayer(prob_surface_split, prob_surface_split_k)
      
      # add to combined surface matrix
      prob_surface_mat <- prob_surface_mat + prob_surface_split_mat/K[i]
    }
    
    # make combined raster
    prob_surface <- setValues(raster_empty, prob_surface_mat)
    
    # produce geoprofile rasters
    geoprofile_split <- raster()
    geoprofile_mat <- 0
    for (k in 1:K[i]) {
      
      # make geoprofile matrix from probability surface
      geoprofile_split_mat <- rank(values(prob_surface_split[[k]]), ties.method = "first")
      geoprofile_split_mat <- 100 * (1 - geoprofile_split_mat/max(geoprofile_split_mat, na.rm = TRUE))
      
      # add raster layer
      geoprofile_split_k <- setValues(raster_empty, geoprofile_split_mat)
      geoprofile_split <- addLayer(geoprofile_split, geoprofile_split_k)
    }
    
    # make combined raster
    geoprofile_mat <- rank(values(prob_surface), ties.method = "first")
    geoprofile_mat <- 100 * (1 - geoprofile_mat/max(geoprofile_mat, na.rm = TRUE))
    geoprofile <- setValues(raster_empty, geoprofile_mat)
    
    # get whether rungs have converged
    converged <- output_raw[[i]]$rung_converged
    if (all_converged && any(!converged)) {
      all_converged <- FALSE
    }
    
    # ---------- ESS ----------
    
    # get ESS
    ESS <- effectiveSize(loglike_sampling)
    ESS[ESS == 0] <- samples # if no variation then assume zero autocorrelation
    ESS[ESS > samples] <- samples # ESS cannot exceed actual number of samples taken
    names(ESS) <- rung_names
    
    # ---------- model comparison statistics ----------
    mu <- mean(loglike_sampling[,ncol(loglike_sampling)])
    sigma_sq <- var(loglike_sampling[,ncol(loglike_sampling)])
    DIC_gelman <- mu + sigma_sq/2
    
    # ---------- acceptance rates ----------
    
    # process acceptance rates
    source_accept <- output_raw[[i]]$source_accept/samples
    names(source_accept) <- group_names
    
    sigma_accept <- output_raw[[i]]$sigma_accept/samples
    names(sigma_accept) <- group_names
    
    #coupling_accept <- output_raw[[i]]$coupling_accept/samples
    
    # ---------- save arguments ----------
    
    output_args <- list(precision_lat = precision_lat,
                        precision_lon = precision_lon,
                        burnin = burnin,
                        samples = samples,
                        auto_converge = auto_converge,
                        converge_test = converge_test,
                        solve_label_switching_on = solve_label_switching_on,
                        pb_markdown = pb_markdown,
                        silent = silent)
    
    # ---------- save results ----------
    
    # add to project
    project$output$single_set[[s]]$single_K[[K[i]]] <- list()
    
    project$output$single_set[[s]]$single_K[[K[i]]]$summary <- list(loglike_intervals = loglike_intervals,
                                                                    prob_surface_split = prob_surface_split,
                                                                    prob_surface = prob_surface,
                                                                    geoprofile_split = geoprofile_split,
                                                                    geoprofile = geoprofile,
                                                                    qmatrix = qmatrix,
                                                                    sigma_intervals = sigma_intervals,
                                                                    expected_popsize_intervals = expected_popsize_intervals,
                                                                    ESS = ESS,
                                                                    DIC_gelman = DIC_gelman,
                                                                    converged = converged,
                                                                    source_accept = source_accept,
                                                                    sigma_accept = sigma_accept)
    
    if (store_raw) {
      project$output$single_set[[s]]$single_K[[K[i]]]$raw <- list(loglike_burnin = loglike_burnin,
                                                                  loglike_sampling = loglike_sampling,
                                                                  source_lon = full_source_lon,
                                                                  source_lat = full_source_lat,
                                                                  sigma = full_sigma,
                                                                  expected_popsize = full_expected_popsize)
    }
    
    project$output$single_set[[s]]$single_K[[K[i]]]$function_call <- list(args = output_args,
                                                                          call = match.call())
    
  } # end loop over K
  
  # name output over K
  K_all <- length(project$output$single_set[[s]]$single_K)
  names(project$output$single_set[[s]]$single_K) <- paste0("K", 1:K_all)
  
  # ---------- tidy up and end ----------
  
  # reorder qmatrices
  project <- align_qmatrix(project)
  
  # run ring-search prior to MCMC
  ringsearch <- ring_search(project,
                            precision_lon = precision_lon,
                            precision_lat = precision_lat,
                            min_lon_model = min_lon_model,
                            max_lon_model = max_lon_model,
                            min_lat_model = min_lat_model,
                            max_lat_model = max_lat_model)
  project$output$single_set[[s]]$all_K$ringsearch <- ringsearch
  
  # get DIC over all K
  DIC_gelman <- mapply(function(x) {
                        ret <- x$summary$DIC_gelman
                        if (is.null(ret)) {
                          return(NA)
                        } else {
                          return(ret)
                        }
                      }, project$output$single_set[[s]]$single_K)
  DIC_gelman <- as.vector(unlist(DIC_gelman))
  project$output$single_set[[s]]$all_K$DIC_gelman <- data.frame(K = 1:length(DIC_gelman), DIC_gelman = DIC_gelman)
  
  # end timer
  tdiff <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  if (tdiff < 60) {
    message(sprintf("Total run-time: %s seconds", round(tdiff, 2)))
  } else {
    message(sprintf("Total run-time: %s minutes", round(tdiff/60, 2)))
  }
  
  # warning if any rungs in any MCMCs did not converge
  if (!all_converged && !silent) {
    message("\n**WARNING** at least one MCMC run did not converge\n")
  }
  
  # return invisibly
  invisible(project)
}

#------------------------------------------------
# align qmatrices over all K
#' @noRd
align_qmatrix <- function(project) {
  
  # get active set
  s <- project$active_set
  
  # extract objects of interest
  x <- project$output$single_set[[s]]$single_K
  
  # find values with output
  null_output <- mapply(function(y) {is.null(y$summary$qmatrix)}, x)
  w <- which(!null_output)
  
  # set template to first qmatrix
  template_qmatrix <- x[[w[1]]]$summary$qmatrix
  n <- nrow(template_qmatrix)
  c <- ncol(template_qmatrix)
  positive_sentinels <- which(!is.na(template_qmatrix[,1]))
  
  # loop through output
  best_perm <- NULL
  for (i in w) {
    
    # expand template
    qmatrix <- unclass(x[[i]]$summary$qmatrix)
    template_qmatrix <- cbind(template_qmatrix, matrix(0, n, i-c))
    
    # calculate cost matrix
    cost_mat <- matrix(0,i,i)
    for (k1 in 1:i) {
      for (k2 in 1:i) {
        cost_mat[k1,k2] <- sum(qmatrix[positive_sentinels,k1] * (log(qmatrix[positive_sentinels,k1]+1e-100) - log(template_qmatrix[positive_sentinels,k2]+1e-100)))
      }
    }
    
    # get lowest cost permutation
    best_perm <- call_hungarian(cost_mat)$best_matching + 1
    best_perm_order <- order(best_perm)
    
    # reorder qmatrix
    group_names <- paste0("group", 1:ncol(qmatrix))
    qmatrix <- qmatrix[, best_perm_order, drop = FALSE]
    colnames(qmatrix) <- group_names
    
    # reorder raw output
    if (!is.null(x[[i]]$raw)) {
      
      # reorder source_lon
      source_lon <- x[[i]]$raw$source_lon[, best_perm_order, drop = FALSE]
      names(source_lon) <- group_names
      project$output$single_set[[s]]$single_K[[i]]$raw$source_lon <- source_lon
      
      # reorder source_lat
      source_lat <- x[[i]]$raw$source_lat[, best_perm_order, drop = FALSE]
      names(source_lat) <- group_names
      project$output$single_set[[s]]$single_K[[i]]$raw$source_lat <- source_lat
      
      # reorder sigma
      sigma <- x[[i]]$raw$sigma
      if (ncol(sigma) > 1) {
        sigma <- sigma[, best_perm_order, drop = FALSE]
        names(sigma) <- group_names
        project$output$single_set[[s]]$single_K[[i]]$raw$sigma <- sigma
      }
    }
    
    # reorder prob_surface_split
    prob_surface_split <- x[[i]]$summary$prob_surface_split
    layer_names <- names(prob_surface_split)
    prob_surface_split <- prob_surface_split[[best_perm_order]]
    names(prob_surface_split) <- layer_names
    project$output$single_set[[s]]$single_K[[i]]$summary$prob_surface_split <- prob_surface_split
    
    # reorder geoprofile_split
    geoprofile_split <- x[[i]]$summary$geoprofile_split
    layer_names <- names(geoprofile_split)
    geoprofile_split <- geoprofile_split[[best_perm_order]]
    names(geoprofile_split) <- layer_names
    project$output$single_set[[s]]$single_K[[i]]$summary$geoprofile_split <- geoprofile_split
    
    # reorder sigma_intervals
    sigma_intervals <- x[[i]]$summary$sigma_intervals[best_perm_order,,drop = FALSE]
    rownames(sigma_intervals) <- group_names
    project$output$single_set[[s]]$single_K[[i]]$summary$sigma_intervals <- sigma_intervals
    
    # reorder source_accept
    source_accept <- x[[i]]$summary$source_accept[best_perm_order]
    names(source_accept) <- group_names
    project$output$single_set[[s]]$single_K[[i]]$summary$source_accept <- source_accept
    
    # reorder sigma_accept
    sigma_accept <- x[[i]]$summary$sigma_accept[best_perm_order]
    names(sigma_accept) <- group_names
    project$output$single_set[[s]]$single_K[[i]]$summary$sigma_accept <- sigma_accept
    
    # qmatrix becomes template for next level up
    template_qmatrix <- qmatrix
    
    # store result
    class(qmatrix) <- "rgeoprofile_qmatrix"
    project$output$single_set[[s]]$single_K[[i]]$summary$qmatrix <- qmatrix
  }
  
  # return modified project
  return(project)
}

#------------------------------------------------
# ring-search
#' @noRd
ring_search <- function(project,
                        precision_lon,
                        precision_lat,
                        min_lon_model,
                        max_lon_model,
                        min_lat_model,
                        max_lat_model) {
  
  # extract sentinel locations with at least one observation
  data <- subset(project$data, counts > 0)
  sentinel_lon <- data$longitude
  sentinel_lat <- data$latitude
  
  # get breaks, midpoints, and coordinates of all points in grid
  breaks_lon <- seq(min_lon_model, max_lon_model, precision_lon)
  breaks_lat <- seq(min_lat_model, max_lat_model, precision_lat)
  midpoints_lon <- breaks_lon[-1] - precision_lon/2
  midpoints_lat <- breaks_lat[-1] - precision_lat/2
  coords <- expand.grid(midpoints_lon, midpoints_lat)
  
  # get distance between all cells and sentinels
  d <- apply(cbind(sentinel_lon, sentinel_lat), 1,
             function(x) lonlat_to_bearing(coords[,1], coords[,2], x[1], x[2])$gc_dist)
  
  # get minimum distance to any sentinel
  d_min <- apply(d, 1, min)
  
  # convert min distances to hitscore percentages
  hs <- rank(d_min)/length(d_min) * 100
  
  # get into raster format
  ret <- raster(xmn = min_lon_model, xmx = max_lon_model,
                ymn = min_lat_model, ymx = max_lat_model,
                nrow = length(midpoints_lat), ncol = length(midpoints_lon))
  ret <- setValues(ret, hs)
  ret <- flip(ret, 2)
  
  return(ret)
}