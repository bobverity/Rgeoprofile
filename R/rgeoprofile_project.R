
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
  
  message("TODO/n")
  
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
