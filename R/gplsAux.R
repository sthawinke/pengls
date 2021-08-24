#' Extract coefficients from a gpls model
#' @param object A gpls object
#' @param ... further arguments, currently ignored
#'
#' @method coef gpls
#' @export
coef.gpls = function(object, ...) {
    coef(object$glmnet)
}
#' Make predictions from a gpls model
#' @param object A gpls object
#' @param ... further arguments, currently ignored
#' @export
#' @method predict gpls
predict.gpls = function(object, ...) {
    predict(object$glmnet, newx = object$data)
}
#' Print a summary of a gpls model
#' @param x A gpls object
#' @param ... further arguments, currently ignored
#' @export
#' @method print gpls
print.gpls = function(x, ...) {
    cat("GPLS model with", sum(coef(x)!=0)-1, "non-zero coefficients,\n and with", corStruct)
}
