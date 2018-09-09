
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
#' @import fftwtools
#' @import ggplot2
#' @import ggmap
#' @import RColorBrewer
#' @import rgdal
#' @importFrom raster raster extent extent<- rasterize projectRaster distance
#' @import viridis
#' @importFrom grDevices colorRampPalette grey
#' @import graphics
#' @import stats
#' @import utils
NULL

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
#' @param source_min_lon the minimum possible longitude of source locations
#' @param source_max_lon the minimum possible longitude of source locations
#' @param source_min_lat the minimum possible latitude of source locations
#' @param source_max_lat the minimum possible latitude of source locations
#' @param sigma_model set as "single" to assume the same dispersal distance for
#'   all sources, or "independent" to assume an independently drawn dispersal
#'   distance for each source
#' @param sigma_prior_mean the prior mean of the parameter sigma (km)
#' @param sigma_prior_sd the prior standard deviation of the parameter sigma 
#'   (km). Set to zero to use a fixed value for sigma
#' @param expected_popsize_prior_mean the prior mean of the expected total
#'   population size
#' @param expected_popsize_prior_sd the prior standard deviation of the expected
#'   total population size
#' 
#' @export

new_set <- function(project,
                    name = "(no name)",
                    sentinel_radius = 0.2,
                    source_min_lon = -10,
                    source_max_lon = 10,
                    source_min_lat = -10,
                    source_max_lat = 10,
                    sigma_model = "single",
                    sigma_prior_mean = 1,
                    sigma_prior_sd = 1,
                    expected_popsize_prior_mean = 100,
                    expected_popsize_prior_sd = 10) {
  
  # check inputs
  assert_custom_class(project, "rgeoprofile_project")
  assert_single_pos(sentinel_radius, zero_allowed = FALSE)
  assert_single_numeric(source_min_lon)
  assert_single_numeric(source_max_lon)
  assert_gr(source_max_lon, source_min_lon)
  assert_single_numeric(source_min_lat)
  assert_single_numeric(source_max_lat)
  assert_gr(source_max_lat, source_min_lat)
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
                                      source_min_lon = source_min_lon,
                                      source_max_lon = source_max_lon,
                                      source_min_lat = source_min_lat,
                                      source_max_lat = source_max_lat,
                                      sigma_model = sigma_model,
                                      sigma_prior_mean = sigma_prior_mean,
                                      sigma_prior_sd = sigma_prior_sd,
                                      expected_popsize_prior_mean = expected_popsize_prior_mean,
                                      expected_popsize_prior_sd = expected_popsize_prior_sd)
  
  names(project$parameter_sets)[s] <- paste0("set", s)
  
  project$output$single_set[[s]] <- list(single_K = list())
  
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
#'   longer be correct, although evidence estimates will be unaffected.
#' @param cluster option to pass in a cluster environment (see package 
#'   "parallel")
#' @param pb_markdown whether to run progress bars in markdown mode, in which 
#'   case they are updated once at the end to avoid large amounts of output.
#' @param silent whether to suppress all console output
#' 
#' @export

run_mcmc <- function(project, K = 3, precision_lon = 1e-3, precision_lat = 1e-3, burnin = 1e2, samples = 1e3, auto_converge = TRUE, converge_test = ceiling(burnin/10), solve_label_switching_on = TRUE, cluster = NULL, pb_markdown = FALSE, silent = !is.null(cluster)) {
  
  # start timer
  t0 <- Sys.time()
  
  # check inputs
  assert_custom_class(project, "rgeoprofile_project")
  assert_single_pos_int(K, zero_allowed = FALSE)
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
  if (s==0) {
    stop("no active parameter set")
  }
  
  # get useful quantities
  n_trap <- nrow(project$data)
  
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
  
  return(output_raw)
  
  #------------------------
  
  # begin processing results
  if (!silent) {
    cat("Processing results\n")
  }
  
  # loop through K
  ret <- list()
  for (i in 1:length(K)) {
    
    # create name lists
    #ind_names <- paste0("ind", 1:n)
    #locus_names <- paste0("locus", 1:L)
    #deme_names <- paste0("deme", 1:K[i])
    #rung_names <- paste0("rung", 1:rungs)
    
    # ---------- raw mcmc results ----------
    
    # get loglikelihood in coda::mcmc format
    loglike_burnin <- mapply(function(x){mcmc(x)}, output_raw[[i]]$loglike_burnin)
    loglike_sampling <- mcmc(t(rcpp_to_mat(output_raw[[i]]$loglike_sampling)))
    
    # alpha
    alpha <- NULL
    if (admix_on) {
      alpha <- mcmc(output_raw[[i]]$alpha_store)
    }
    
    # ---------- summary results ----------
    
    # get quantiles over sampling loglikelihoods
    loglike_quantiles <- t(apply(loglike_sampling, 2, quantile_95))
    rownames(loglike_quantiles) <- rung_names
    class(loglike_quantiles) <- "maverick_loglike_quantiles"
    
    # process qmatrix_ind
    qmatrix_ind <- rcpp_to_mat(output_raw[[i]]$qmatrix_ind)
    colnames(qmatrix_ind) <- deme_names
    rownames(qmatrix_ind) <- ind_names
    class(qmatrix_ind) <- "maverick_qmatrix_ind"
    
    # ---------- ESS ----------
    
    # get ESS
    ESS <- effectiveSize(loglike_sampling)
    ESS[ESS == 0] <- samples # if no variation then assume zero autocorrelation
    ESS[ESS > samples] <- samples # ESS cannot exceed actual number of samples taken
    names(ESS) <- rung_names
    
    # ---------- acceptance rates ----------
    
    # process acceptance rates
    coupling_accept <- output_raw[[i]]$coupling_accept/samples
    
    # ---------- save arguments ----------
    
    output_args <- list(burnin = burnin,
                        samples = samples,
                        rungs = rungs,
                        GTI_pow = GTI_pow,
                        auto_converge = auto_converge,
                        converge_test = converge_test,
                        solve_label_switching_on = solve_label_switching_on,
                        coupling_on = coupling_on,
                        pb_markdown = pb_markdown,
                        silent = silent)
    
    # ---------- save results ----------
    
    # add to project
    project$output$single_set[[s]]$single_K[[K[i]]] <- list()
    
    project$output$single_set[[s]]$single_K[[K[i]]]$summary <- list(qmatrix_ind = qmatrix_ind,
                                                                    loglike_quantiles = loglike_quantiles,
                                                                    ESS = ESS,
                                                                    GTI_path = GTI_path,
                                                                    GTI_logevidence = GTI_logevidence)
    
    project$output$single_set[[s]]$single_K[[K[i]]]$raw <- list(loglike_burnin = loglike_burnin,
                                                                loglike_sampling = loglike_sampling,
                                                                alpha = alpha,
                                                                coupling_accept = coupling_accept)
    
    project$output$single_set[[s]]$single_K[[K[i]]]$function_call <- list(args = output_args,
                                                                          call = match.call())
    
  } # end loop over K
  
  # name output over K
  K_all <- length(project$output$single_set[[s]]$single_K)
  names(project$output$single_set[[s]]$single_K) <- paste0("K", 1:K_all)
  
  # ---------- tidy up and end ----------
  
  # reorder qmatrices
  #project <- align_qmatrix(project)
  
  # end timer
  tdiff <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  if (tdiff<60) {
    message(sprintf("Total run-time: %s seconds", round(tdiff, 2)))
  } else {
    message(sprintf("Total run-time: %s minutes", round(tdiff/60, 2)))
  }
  
  # return invisibly
  invisible(project)
}
