#'Iterative estimation of penalised generalised least squares
#'
#' @param data A data matrix
#' @param glsSt a covariance structure, as supplied to nlme::gls as "correlation"
#' @param optControl control arguments, passed onto nlme::gls' control argument
#' @param xNames names of the regressors in data
#' @param outVar name of the outcome variable in data
#' @param coords Names of the spatial covariates, used by glsSt
#' @param corMat a starting value for th correlation matrix. Taken to be a diagonal matrix if missing
#' @param Coef a starting value for coefficients
#' @param maxIter maximum number of iterations between glmnet and gls
#' @param tol A convergence tolerance
#' @param verbose a boolean, should output be printed?
#' @param ... passed onto glmnet::glmnet
#' @return
#' @import glmnet
pgls = function(data, glsSt, xNames, outVar, coords = c("x", "y"), optControl = list(), corMat,
                     maxIter = 3e1, tol = 5e-2, verbose = FALSE, ...){
   iter = 1L; conv = FALSE
   if(missing(corMat))
        corMat = diag(nrow(data)) #Starting values for correlation matrix
    preds = mcA = data[[outVar]] - mean(data[[outVar]]) #Starting values for predictions
    xY = as.matrix(data[,c(outVar, xNames)]) #The outcome and design matrix
    desMat = cbind(1, xY[, -1]) #The design matrix
    while(iter <= maxIter && !conv){
        oldPred = preds #Store old predictions
        tmpDat = corMat %*% xY #Pre-multiply data by correlation matrix
        glmnetFit = glmnet(x = tmpDat[,-1], y = tmpDat[,1], intercept = TRUE, ...) #Fit glmnet
        preds = drop(meanCenterVec(desMat %*% coef(glmnetFit))) #Make predictions and center
        margCorMat = getCorMat(kernelData = cbind("a" = drop(mcA - preds), data[, c("x", "y")]),
                               control = optControl, glsSt = glsSt)#Find the correlation matrix
        corMat = margCorMat$corMat;Coef = margCorMat$Coef
        resSqt = sqrt(mean(((preds-oldPred))^2)) #The mean squared change in predictions
        conv = resSqt < tol #Check for convergence
        iter = iter + 1L
    }
    if(!conv && verbose) {
        warning("No convergence achieved in PGLS!\n", immediate. = TRUE)
    }
    glmnetFit$corMat = corMat
    return(glmnetFit)
}
