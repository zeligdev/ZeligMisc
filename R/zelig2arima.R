#' Interface between the Zelig Model arima and 
#' the Pre-existing Model-fitting Method
#' @param formula a formula
#' @param ... additonal parameters
#' @param data a data.frame 
#' @return a list specifying '.function'
#' @export
zelig2arima <- function (formula, ..., data) {



  list(
       # "arima.wrap" 
       .function = "arima.wrap",
       formula = formula,
       data = data
       )
}
