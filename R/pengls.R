#'Iterative estimation of penalised generalised least squares
#'
#' @param data A data matrix or data frame
#' @param glsSt a covariance structure, as supplied to nlme::gls as "correlation"
#' @param optControl control arguments, passed onto nlme::gls' control argument
#' @param xNames names of the regressors in data
#' @param outVar name of the outcome variable in data
#' @param corMat a starting value for the correlation matrix. Taken to be a diagonal matrix if missing
#' @param lambda The penalty value for glmnet. If missing, the optimal value of vanilla glmnet without autocorrelation component is used
#' @param foldid An optional vector defining the fold
#' @param exclude indices of predictors to be excluded from intercept + xNames
#' @param maxIter maximum number of iterations between glmnet and gls
#' @param tol A convergence tolerance
#' @param verbose a boolean, should output be printed?
#' @param scale,center booleans, should regressors be scaled to zero mean and variance 1? Defaults to TRUE
#' @param nfolds an integer, the number of folds used in cv.glmnet to find lambda
#' @param penalty.factor passed onto glmnet:glmnet. The first entry is zero by default for the intercept, which is not shrunk
#' @param ... passed onto glmnet::glmnet
#' @return A list with components
#' \item{glmnet}{The glmnet fit, which can be manipulated as such}
#' \item{gls}{A list with info on the estimated correlation matrix}
#' \item{iter}{The iterations needed}
#' \item{conv}{A boolean, indicating whether the iteration between mean model and covariance estimation converged}
#' \item{xNames,data,glsSt,outVar}{As provided}
#' \item{lambda}{The lambda penalty paraneter used}
#' @import glmnet nlme
#' @importFrom stats coef
#' @seealso cv.pengls
#' @export
#' @examples
#' ### Example 1: spatial data
#' # Define the dimensions of the data
#' library(nlme)
#' n <- 50 #Sample size
#' p <- 100 #Number of features
#' g <- 10 #Size of the grid
#' #Generate grid
#' Grid <- expand.grid("x" = seq_len(g), "y" = seq_len(g))
#' # Sample points from grid without replacement
#' GridSample <- Grid[sample(nrow(Grid), n, replace = FALSE),]
#' #Generate outcome and regressors
#' b <- matrix(rnorm(p*n), n , p)
#' a <- rnorm(n, mean = b %*% rbinom(p, size = 1, p = 0.2)) #20% signal
#' #Compile to a matrix
#' df <- data.frame("a" = a, "b" = b, GridSample)
#' # Define the correlation structure (see ?nlme::gls), with initial nugget 0.5 and range 5
#' corStruct <- corGaus(form = ~ x + y, nugget = TRUE, value = c("range" = 5, "nugget" = 0.5))
#' #Fit the pengls model, for simplicity for a simple lambda
#' penglsFit <- pengls(data = df, outVar = "a", xNames = grep(names(df), pattern = "b", value =TRUE),
#' glsSt = corStruct, nfolds = 5)
#'
#' ### Example 2: timecourse data
#' dfTime <- data.frame("a" = a, "b" = b, "t" = seq_len(n))
#' dfTime$a[-1] = dfTime$a[-n]*0.25 #Some temporal signal
#' corStructTime <- corAR1(form = ~ t, value = 0.5)
#' penglsFitTime <- pengls(data = dfTime, outVar = "a",
#' xNames = grep(names(dfTime), pattern = "b", value =TRUE),
#' glsSt = corStructTime, nfolds = 5)
pengls = function(data, glsSt, xNames, outVar, corMat, lambda, foldid, exclude = NULL,
                maxIter = 3e1, tol = 5e-2, verbose = FALSE, scale = FALSE, center = FALSE,
                optControl = lmeControl(opt = "optim", maxIter = 5e2, msVerbose = verbose,
                                        msMaxIter = 5e2, niterEM = 1e3,
                                        msMaxEval=1e3), nfolds = 10, penalty.factor = c(0, rep(1, length(xNames))), ...){
   coords <- {
      foo = trimws(strsplit(split = "\\+", as.character(attr(glsSt, "formula"))[2])[[1]]) #Extract the coordinates from the formula
      foo[!foo %in% c("+", " ")]
   }
   glsSt0 = glsStruct(glsSt) #Initialize correlation object
   if(missing(corMat) )
        corMat <- diag(nrow(data)) #Starting values for correlation matrix
   if(any(c(scale, center)))
      data[, xNames] = scale(data[, xNames], scale = scale, center = center) #Center and scale
   if(missing(lambda)){
       if(verbose) cat("Fitting naive model...\n")
       naieveFit <- cv.glmnet(x = as.matrix(data[,xNames]), y = data[,outVar], exclude = exclude-1,
                              nfolds = nfolds, penalty.factor = penalty.factor[-1], ...) #Exclude intercept heres
       lambda <- naieveFit$lambda.1se
   }
    preds <- mcA <- data[[outVar]] - mean(data[[outVar]]) #Starting values for predictions
    xY = cbind("Intercept" = 1, as.matrix(data[, c(outVar, xNames)]))
    if(verbose) cat("Starting iterations...\n")
    iter = 1L; conv = FALSE
    while(iter <= maxIter && !conv){
        if(verbose) cat("Iteration", iter, "\n")
        oldPred <- preds #Store old predictions
        tmpDat <- corMat %*% xY #Pre-multiply data by correlation matrix
        glmnetFit <- glmnet(x = tmpDat[,c("Intercept", xNames)], y = tmpDat[,outVar],
                            intercept = FALSE, penalty.factor = penalty.factor, exclude = exclude,
                            lambda = lambda, ...) #Fit glmnet, do not shrink intercept
        preds <- as.vector(xY[, -2] %*% coef(glmnetFit)[-1]) #Make predictions
        margCorMat <- getCorMat(data = cbind("a" = mcA - preds, data[, coords, drop = FALSE]), outVar = "a", control = optControl, glsSt = glsSt0)#Find the correlation matrix
        corMat <- margCorMat$corMat;Coef = margCorMat$Coef
        resSqt <- sqrt(mean(((preds-oldPred))^2)) #The mean squared change in predictions
        conv <- resSqt < tol #Check for convergence
        iter <- iter + 1L
    }
    if(!conv && verbose) {
        warning("No convergence achieved in pengls!\n", immediate. = TRUE)
    }
    out <- list("glmnet" = glmnetFit, "gls" = margCorMat, "data" = data,
               "xNames" = xNames, "outVar" = outVar, "glsSt" = glsSt,
               "lambda" = lambda, "iter" = iter, "conv" = conv)
    class(out) <- "pengls"
    return(out)
}
