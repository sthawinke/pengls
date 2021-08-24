#' Extract coefficients from a gpls model
#' @param object A gpls object
#' @param ...
#'
#' @method coef gpls
#' @export
gpls.coef = function(object, ...) {
    coef(object$glmnet)
}
