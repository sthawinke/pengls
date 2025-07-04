#' Extract coefficients from a pengls model
#' @param object A pengls object
#' @param ... further arguments, currently ignored
#'
#' @method coef pengls
#' @export
#' @return The vector of coefficients
coef.pengls = function(object, ...) {
    Coef <- coef(object$glmnet)[-1]
    names(Coef) = c("Intercept", object$xNames)
    return(Coef)
}
#' Extract coefficients from a cv.pengls model
#' @param object A cv.pengls object
#' @param which a character string, for which lambda should coefficients be returned
#' @param ... further arguments, currently ignored
#'
#' @method coef cv.pengls
#' @export
#' @return The vector of coefficients
coef.cv.pengls = function(object, which = "lambda.1se", ...) {
    Coef <- object$coefs[,object$lambda == object[[which]]]
    names(Coef) = c("Intercept", object$xNames)
    return(Coef)
}
#' Make predictions from a pengls model
#' @param object A pengls object
#' @param newx The test data
#' @param ... further arguments, currently ignored
#' @export
#' @importFrom stats predict
#' @method predict pengls
#' @return A vector with predicted values
predict.pengls = function(object, newx, ...) {
    if(missing(newx))
        newx = as.matrix(object$data[, object$xNames])
    predict(object$glmnet, newx = cbind(1, newx), ...)
}
#' Make predictions from a cv.pengls model
#' @param object A cv.pengls object
#' @param ... further arguments, currently ignored
#' @export
#' @importFrom stats predict
#' @method predict cv.pengls
#' @return A vector with predicted values
predict.cv.pengls = function(object, ...) {
    predict.pengls(object$bestFit, ...)
}
#' Print a summary of a pengls model
#' @param x A pengls object
#' @param ... further arguments, currently ignored
#' @export
#' @return Prints output to console
#' @method print pengls
print.pengls = function(x, ...) {
    Class <- class(x$glsSt)[1]
    cat("pengls model with correlation structure:", Class, "\n and", sum(coef(x)!=0)-1, "non-zero coefficients")
}
#' Print a summary of a cv.pengls model
#' @param x A cv.pengls object
#' @param ... further arguments, currently ignored
#' @export
#' @return Prints output to console
#' @method print cv.pengls
print.cv.pengls = function(x, ...) {
    Class  <- class(x$glsSt)[1]
    cat("Cross-validated pengls model with correlation structure:", Class, "\n and", sum(coef(x)!=0)-1, "non-zero coefficients.\n",
        length(unique(x$foldid)), "fold cross-validation yielded an estimated", x$loss, "of", x$cvOpt, ".")
}
