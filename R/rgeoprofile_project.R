
#------------------------------------------------
#' @title Create a new RgeoProfile project
#'
#' @description Create a new RgeoProfile project.
#'
#' @export

rgeoprofile_project <- function() {
  
  # initialise project with default values
  ret <- list(data = NULL,
              parameter_sets = NULL,
              active_set = 0,
              output = list(single_set = list())
              )
  
  # create class and return
  class(ret) <- "rgeoprofile_project"
  return(ret)
}

#------------------------------------------------
#' @title Custom print function for class rgeoprofile_project
#'   
#' @description Custom print function for class rgeoprofile_project, printing a summary 
#'   of the key elements (also equivalent to \code{summary(x)}). To do an 
#'   ordinary \code{print()} of all elements of the project, use the 
#'   \code{print_full()} function.
#'   
#' @param x object of class \code{rgeoprofile_project}
#' @param ... other arguments (ignored)
#'   
#' @export

print.rgeoprofile_project <- function(x, ...) {
  
  # print summary
  summary(x)
  
  # return invisibly
  invisible(x)
}

#------------------------------------------------
#' @title Ordinary print function for class rgeoprofile_project
#'
#' @description Calling \code{print()} on an object of class rgeoprofile_project results
#'   in custom output. This function therefore stands in for the base
#'   \code{print()} function, and is equivalent to running
#'   \code{print(unclass(x))}.
#'
#' @param x object of class \code{rgeoprofile_project}
#' @param ... other arguments passed to \code{print()}
#'
#' @export

print_full <- function(x, ...) {
  
  # check inputs
  assert_custom_class(x, "rgeoprofile_project")
  
  # print un-classed object
  print(unclass(x), ...)
  
  # return invisibly
  invisible(x)
}

#------------------------------------------------
#' @title Print summary for class rgeoprofile_project
#'   
#' @description Overload summary function for class rgeoprofile_project
#'   
#' @param object object of class \code{rgeoprofile_project}
#' @param ... other arguments (ignored)
#'   
#' @export

summary.rgeoprofile_project <- function(object, ...) {
  
  # print data summary
  cat("DATA:\n")
  if (is.null(object$data)) {
    cat("   (none loaded)\n")
  } else {
    
    # extract data properties
    n_sentinel <- nrow(object$data)
    n_positive <- sum(object$data$counts > 0)
    n_obs <- sum(object$data$counts)
    
    cat(sprintf("   sentinel sites = %s (%s positive, %s negative)\n", n_sentinel, n_positive, n_sentinel - n_positive))
    cat(sprintf("   total observations = %s\n", n_obs))
  }
  cat("\n")
  
  # print parameter sets summary
  cat("PARAMETER SETS:\n")
  if (length(object$parameter_sets) == 0) {
    cat("   (none defined)\n")
  } else {
    
    # print names of all sets
    s <- object$active_set
    for (i in 1:length(object$parameter_sets)) {
      
      # star next to active set
      if (i == s) {
        cat(" * ")
      } else {
        cat("   ")
      }
      
      # print name of set
      cat(sprintf("SET%s: %s\n", i, object$parameter_sets[[i]]$name))
    }
    cat("\n")
    
    # get details of active set
    sentinel_radius <- object$parameter_sets[[s]]$sentinel_radius
    spatial_prior <- object$parameter_sets[[s]]$spatial_prior
    min_lon <- round(xmin(spatial_prior), digits = 3)
    max_lon <- round(xmax(spatial_prior), digits = 3)
    min_lat <- round(ymin(spatial_prior), digits = 3)
    max_lat <- round(ymax(spatial_prior), digits = 3)
    cells_lon <- ncol(spatial_prior)
    cells_lat <- nrow(spatial_prior)
    sigma_model <- object$parameter_sets[[s]]$sigma_model
    sigma_prior_mean <- object$parameter_sets[[s]]$sigma_prior_mean
    sigma_prior_sd <- object$parameter_sets[[s]]$sigma_prior_sd
    expected_popsize_prior_mean <- object$parameter_sets[[s]]$expected_popsize_prior_mean
    expected_popsize_prior_sd <- object$parameter_sets[[s]]$expected_popsize_prior_sd
    
    # print details of active set
    cat(sprintf("ACTIVE SET: SET%s\n", s))
    cat(sprintf("   sentinel radius = %s\n", sentinel_radius))
    cat(sprintf("   spatial prior:\n"))
    cat(sprintf("      longitude range = [%s, %s]\n", min_lon, max_lon))
    cat(sprintf("      latitude range = [%s, %s]\n", min_lat, max_lat))
    cat(sprintf("      cells = %s, %s (lon, lat)\n", cells_lon, cells_lat))
    cat(sprintf("   sigma prior:\n"))
    cat(sprintf("      model = %s\n", sigma_model))
    if (sigma_prior_sd == 0) {
      cat(sprintf("      value = %s (exact prior)\n", sigma_prior_mean))
    } else {
      cat(sprintf("      prior mean = %s\n", sigma_prior_mean))
      cat(sprintf("      prior SD = %s\n", sigma_prior_sd))
    }
    cat(sprintf("   expected population size prior:\n"))
    if (expected_popsize_prior_sd == -1) {
      cat(sprintf("      (improper prior)\n"))
    } else if (expected_popsize_prior_sd == 0) {
      cat(sprintf("      value = %s (exact prior)\n", expected_popsize_prior_mean))
    } else {
      cat(sprintf("      prior mean = %s\n", expected_popsize_prior_mean))
      cat(sprintf("      prior SD = %s\n", expected_popsize_prior_sd))
    }
    
  }
  cat("\n")
  
}

#------------------------------------------------
#' @title Determine if object is of class rgeoprofile_project
#'
#' @description Determine if object is of class rgeoprofile_project
#'
#' @param x object of class \code{rgeoprofile_project}
#'
#' @export

is.rgeoprofile_project <- function(x) {
  inherits(x, "rgeoprofile_project")
}
