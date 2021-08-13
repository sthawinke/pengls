#Get the (square root of the inverse of the) correlation matrix
getMargCorMat = function(kernelData, glsSt, Coef = c(coef(glsSt)), control, outVar = "a", corr = FALSE){
    attr(glsSt, "conLin") <- list(Xy = cbind(1,kernelData[, outVar]),
                                  dims = list(N = nrow(kernelData), p = 1, REML = 1),
                                  logLik = 0, sigma = control$sigma, fixedSigma = FALSE)
    glsSt = Initialize(glsSt, kernelData[, c(outVar, "x", "y")]) #Initialize the correlation matrix
    optRes = optim(Coef, function(glsPars) -logLik(glsSt, glsPars), method = control$optimMethod,
          control = list(trace = control$msVerbose, maxit = control$msMaxIter, reltol = control$msTol)) #Find the parameters
    coef(glsSt) <- optRes$par
    corMat = corMatrix(glsSt$corStruct, corr = corr) #The inverse square root glsSt matrix, see ?corMat
    return(list("corMat" = corMat, "Coef" = optRes$par))
}