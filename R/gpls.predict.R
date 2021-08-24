#' Make predictions from a gpls model
#' @param object A gpls object
#' @param ...
#' @export
#' @method predict gpls
gpls.predict = function(object, ...) {
    predict(object$glmnet)
}
