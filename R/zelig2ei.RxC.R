#' Interface between the Zelig Model ei.RxC and 
#' the Pre-existing Model-fitting Method
#' @param formula a formula
#' @param ... additonal parameters
#' @param data a data.frame 
#' @return a list specifying '.function'
#' @export
zelig2ei.RxC <- function (formula, ..., data) {
  list(
       .function = "",
       formula = formula,
       data = data
       )
}
