#' Extract coefficients from a gpls model
#' @param object A gpls object
#' @param ... further arguments, currently ignored
#'
#' @method coef gpls
#' @export
#' @return The vector of coefficients
coef.gpls = function(object, ...) {
    coef(object$glmnet)
}
#' Make predictions from a gpls model
#' @param object A gpls object
#' @param ... further arguments, currently ignored
#' @export
#' @importFrom stats predict
#' @method predict gpls
#' @return A vector with predicted values
predict.gpls = function(object, ...) {
    predict(object$glmnet, newx = as.matrix(object$data[, object$xNames]))
}
#' Print a summary of a gpls model
#' @param x A gpls object
#' @param ... further arguments, currently ignored
#' @export
#' @return Prints output to console
#' @method print gpls
print.gpls = function(x, ...) {
    cat("GPLS model with correlation structure:", class(x$glsSt)[1], "\n and", sum(coef(x)!=0)-1, "non-zero coefficients")
}
