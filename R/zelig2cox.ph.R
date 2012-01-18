#' Interface between the Zelig Model coxph and 
#' the Pre-existing Model-fitting Method
#' @param formula a formula
#' @param ... additonal parameters
#' @param data a data.frame 
#' @return a list specifying '.function'
#' @export
zelig2cox.ph <- function (formula, ..., robust = FALSE, cluster, data) {

  # Error-catching
  if (!is.null(cluster) && !robust)
    stop("If cluster is specified, robust must be TRUE.")

  # If a cluster is present, then we update the formula accordingly
  if (!is.null(cluster)) {
    update.str <- sprintf(". ~ . + cluster(%s)", cluster)
    formula <- update(formula, update.str)
  }

  # Return
  list(
       .function = "coxph",
       .hook = "coxph.hook",
       formula = formula,
       robust = robust,
       ...,
       data = data
       )
}


# Add robust/naive support for Cox Proportional Hazard model
#' @export
coxph.hook <- function (obj, Zall, Call, robust, ...) {
  return(obj)
  if (!is.null(robust)) {
    if (!is.logical(robust))
      stop("invalid input for robust.  Choose either TRUE or FALSE.")

    class(obj) <- if (robust)
      c("coxph.robust", class(obj))

    else
      c("coxph.naive", class(obj))
  }
  else
    class(obj) <- c("coxph.naive", class(obj))

  # Return
  class(obj) <- c(class(obj), "zelig")
  obj
}
