#'Iterative estimation of penalised generalised least squares
#'
#' @param data A data matrix or data frame
#' @param glsSt a covariance structure, as supplied to nlme::gls as "correlation"
#' @param optControl control arguments, passed onto nlme::gls' control argument
#' @param xNames names of the regressors in data
#' @param outVar name of the outcome variable in data

#' @param corMat a starting value for th correlation matrix. Taken to be a diagonal matrix if missing
#' @param lambda The penalty value for glmnet. If missing, the optimal value of vanilla glmnet without autocorrelation component is used
#' @param maxIter maximum number of iterations between glmnet and gls
#' @param tol A convergence tolerance
#' @param verbose a boolean, should output be printed?
#' @param nfolds an integer, the number of folds used in cv.glmnet to find lambda
#' @param ... passed onto glmnet::glmnet
#' @return A list with components
#' \item{corMat}{The square root of the inverse correlation matrix}
#' \item{Coef}{The coefficients of the correlation object}
#' @import glmnet nlme
#' @importFrom stats coef
#' @export
#' @details This function does not provide cross-validation, but rather fits the model
#' for the lambda penalty value provided, or else the optimal lambda value of the vanilla glmnet.
#' Cross-validation needs to be implemented by the user; since the data exhibits autocorrelation
#' this may need to be blocked cross-validation or some other dedicated method.
#' @examples
#' ### Example 1: spatial data
#' # Define the dimensions of the data
#' library(nlme)
#' n = 50 #Sample size
#' p = 100 #Number of features
#' g = 10 #Size of the grid
#' #Generate grid
#' Grid = expand.grid("x" = seq_len(g), "y" = seq_len(g))
#' # Sample points from grid without replacement
#' GridSample = Grid[sample(nrow(Grid), n, replace = FALSE),]
#' #Generate outcome and regressors
#' a = rnorm(n)
#' b = matrix(rnorm(p*n), n , p)
#' #Compile to a matrix
#' df = data.frame("a" = a, "b" = b, GridSample)
#' # Define the correlation structure (see ?nlme::gls), with initial nugget 0.5 and range 5
#' corStruct = corGaus(form = ~ x + y, nugget = TRUE, value = c("range" = 5, "nugget" = 0.5))
#' #Fit the gpls model, for simplicity for a simple lambda
#' gplsFit = gpls(data = df, outVar = "a", xNames = grep(names(df), pattern = "b", value =TRUE),
#' glsSt = corStruct, nfolds = 5)
#'
#' ### Example 2: timecourse data
#' dfTime = data.frame("a" = a, "b" = b, "t" = seq_len(50))
#' corStructTime = corAR1(form = ~ t, value = 0.5)
#' gplsFit = gpls(data = dfTime, outVar = "a",
#' xNames = grep(names(dfTime), pattern = "b", value =TRUE),
#' glsSt = corStructTime, nfolds = 5)
gpls = function(data, glsSt, xNames, outVar, corMat, lambda,
                     maxIter = 3e1, tol = 5e-2, verbose = FALSE,
                optControl = lmeControl(opt = "optim", maxIter = 5e2, msVerbose = verbose,
                                        msMaxIter = 5e2, niterEM = 1e3,
                                        msMaxEval=1e3), nfolds = 10,  ...){
   iter = 1L; conv = FALSE
   if(missing(corMat))
        corMat = diag(nrow(data)) #Starting values for correlation matrix
   if(missing(lambda)){
       if(verbose) cat("Fitting naieve model...\n")
       naieveFit = cv.glmnet(x = as.matrix(data[,xNames]), y = data[,outVar],
                             nfolds = nfolds, ...)
       lambda = naieveFit$lambda.1se
   }
    preds = mcA = data[[outVar]] - mean(data[[outVar]]) #Starting values for predictions
    xY = as.matrix(data[,c(outVar, xNames)]) #The outcome and design matrix
    desMat = cbind(1, xY[, -1]) #The design matrix
    if(verbose) cat("Starting iterations...\n")
    while(iter <= maxIter && !conv){
        if(verbose) cat("Iteration", iter, "\n")
        oldPred = preds #Store old predictions
        tmpDat = corMat %*% xY #Pre-multiply data by correlation matrix
        glmnetFit = glmnet(x = tmpDat[,xNames], y = tmpDat[,outVar], lambda = lambda, ...) #Fit glmnet
        preds = as.vector(desMat %*% coef(glmnetFit)) #Make predictions and center
        margCorMat = getCorMat(data = cbind("a" = mcA - preds, data), outVar = outVar,
                               control = optControl, glsSt = glsSt)#Find the correlation matrix
        corMat = margCorMat$corMat;Coef = margCorMat$Coef
        resSqt = sqrt(mean(((preds-oldPred))^2)) #The mean squared change in predictions
        conv = resSqt < tol #Check for convergence
        iter = iter + 1L
    }
    if(!conv && verbose) {
        warning("No convergence achieved in gpls!\n", immediate. = TRUE)
    }
    out = list("glmnet" = glmnetFit, "gls" = margCorMat, "data" = data,
               "xNames" = xNames, "outVar" = outVar, "glsSt" = glsSt)
    class(out) = "gpls"
    return(out)
}
