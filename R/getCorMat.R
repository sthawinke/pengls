#' Get the (square root of the inverse of the) correlation matrix
#'
#' @param data The data frame
#' @param glsSt The correlation object for gls
#' @param Coef optional vector of coefficients to glsSt
#' @param control the list of control arguments for gls
#' @param outVar the name of the outcome variable
#'
#' @return A list with components
#' \item{corMat}{The square root of the inverse correlation matrix}
#' \item{Coef}{The coefficients of the correlation object}
#' @importFrom nlme Initialize corMatrix
#' @importFrom stats coef optim logLik
getCorMat = function(data, glsSt, Coef = c(coef(glsSt)), control, outVar){
    attr(glsSt, "conLin") <- list(Xy = cbind(1,data[, outVar]),
                                  dims = list(N = nrow(data), p = 1, REML = 1),
                                  logLik = 0, sigma = control$sigma, fixedSigma = FALSE)
    glsSt <- Initialize(glsSt, data) #Initialize the correlation matrix
    optRes <- optim(Coef, function(glsPars) -logLik(glsSt, glsPars), method = control$optimMethod,
          control = list(maxit = control$msMaxIter, reltol = control$msTol)) #Find the parameters
    coef(glsSt) <- optRes$par
    corMat <- corMatrix(glsSt, corr = FALSE) #The inverse square root glsSt matrix, see ?corMat
    return(list("corMat" = corMat, "Coef" = optRes$par))
}
