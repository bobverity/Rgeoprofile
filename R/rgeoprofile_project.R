
#------------------------------------------------
#' @title Create a new RgeoProfile project
#'
#' @description Create a new RgeoProfile project.
#'
#' @export

rgeoprofile_project <- function() {
  
  # initialise project with default values
  ret <- list(data = NULL,
              data_processed = NULL,
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
    name <- object$data_processed$name
    ploidy <- object$data_processed$ploidy
    ploidy_min <- min(ploidy)
    ploidy_max <- max(ploidy)
    n <- length(ploidy)
    loci <- length(object$data_processed$Jl)
    pop <- object$data_processed$pop
    
    if (!is.null(name)) {
      cat(sprintf("   '%s'\n", name))
    }
    cat(sprintf("   individuals = %s\n", n))
    cat(sprintf("   loci = %s\n", loci))
    if (ploidy_min==ploidy_max) {
      cat(sprintf("   ploidy = %s\n", ploidy_min))
    } else {
      cat(sprintf("   ploidy range = %s to %s\n", ploidy_min, ploidy_max))
    }
    if (is.null(pop)) {
      cat("   pops = (none defined)\n")
    } else {
      cat(sprintf("   pops = %s\n", length(unique(pop))))
    }
    
    n1 <- sum(object$data_processed$dat==0)
    n2 <- length(object$data_processed$dat)
    cat(sprintf("   missing data = %s of %s gene copies (%s%%)\n", n1, n2, round(n1/n2*100)))
  }
  cat("\n")
  
  # print parameter sets summary
  cat("PARAMETER SETS:\n")
  if (length(object$parameter_sets)==0) {
    cat("   (none defined)\n")
  } else {
    # print names of all sets
    s <- object$active_set
    for (i in 1:length(object$parameter_sets)) {
      
      # star next to active set
      if (i==s) {
        cat(" * ")
      } else {
        cat("   ")
      }
      
      # print name of set
      cat(sprintf("SET%s: %s\n", i, object$parameter_sets[[i]]$name))
    }
    cat("\n")
    
    # print details of active set
    name <- object$parameter_sets[[s]]$name
    admix_on <- object$parameter_sets[[s]]$admix_on
    estimate_alpha <- object$parameter_sets[[s]]$estimate_alpha
    alpha <- object$parameter_sets[[s]]$alpha
    lambda <- object$parameter_sets[[s]]$lambda
    
    cat(sprintf("ACTIVE SET: SET%s\n", s))
    if (admix_on) {
      cat("   model = admixture\n")
      cat(sprintf("   estimate alpha = %s\n", estimate_alpha))
      if (!estimate_alpha) {
        cat(sprintf("   alpha = %s\n", alpha))
      }
    } else {
      cat("   model = no-admixture\n")
    }
    cat(sprintf("   lambda = %s\n", lambda))
    
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
