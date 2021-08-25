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
#' Extract coefficients from a cv.gpls model
#' @param object A cv.gpls object
#' @param which a character string, for which lambda shoudl coefficients be returned
#' @param ... further arguments, currently ignored
#'
#' @method coef cv.gpls
#' @export
#' @return The vector of coefficients
coef.cv.gpls = function(object, which = "lambda.1se", ...) {
    object$coefs[,object$lambda == object[[which]]]
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
    Class = class(x$glsSt)[1]
    cat("GPLS model with correlation structure:", Class, "\n and", sum(coef(x)!=0)-1, "non-zero coefficients")
}
#' Print a summary of a cv.gpls model
#' @param x A cv.gpls object
#' @param ... further arguments, currently ignored
#' @export
#' @return Prints output to console
#' @method print cv.gpls
print.cv.gpls = function(x, ...) {
    Class = class(x$glsSt)[1]
    cat("Cross-validated GPLS model with correlation structure:", Class, "\n and", sum(coef(x)!=0)-1, "non-zero coefficients.\n",
        length(unique(x$foldid)), "fold cross-validation yielded an estimated R² of", x$cvOpt, ".")
}
