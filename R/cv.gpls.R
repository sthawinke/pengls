#' Peform cross-validation gpls
#'
#' @inheritParams gpls
#' @param lambdas an optional lambda sequence
#'
#' @return
#' @export
#'
#' @examples
#' #' library(nlme)
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
#' gplsFitCV = cv.gpls(data = df, outVar = "a", xNames = grep(names(df), pattern = "b", value =TRUE),
#' glsSt = corStruct, nfolds = 5)
#' gplsFitCV$lambda.1se #Lambda for 1 standard error rule
#' coef(gplsFitCV)
cv.gpls = function(data, glsSt, xNames, outVar, corMat, foldid, cvType, lambdas, transFun = "identity", transFunArgs = list()){
    if(missing(foldid))
        foldid = makeFolds(nfolds, data, cvType)
    if(missing(lambdas)){
        naieveFit = cv.glmnet(x = as.matrix(data[,xNames]), y = data[,outVar], foldid = foldid, ...)#The naive fit provides a lambda sequence
        lambdas = naieveFit$lambdas
    }
    uniqueFolds  = unique(foldid)
    fits = lapply(uniqueFolds, function(i){
        id = foldid!=i #Leave out test fold
        lapply(lambdas, function(lam){
                gpls(data = data[id,], glsSt = glsSt, xNames = xNames, outVar = outVar, corMat = corMat, lambda = lam)
        })
    })
    #Find R2 on left out folds
    R2sCv = sapply(seq_along(fits), function(i){
        id = foldid==uniqueFolds[i]
        coefs = sapply(fits[[i]], coef)
        preds = apply(cbind(1, as.matrix(data[id,xNames, drop = FALSE])) %*% coefs, 2,
                      function(x) do.call(what = transFun, args = transFunArgs, x))
        testDat = do.call(what = transFun, args = transFunArgs, data[id,outVar])
        1-colMeans((preds-testDat)^2)/varN(testDat)
    })
    #Use lambda within 1se from minimum error
    cvEsts = rowMeans(R2sCv)
    maxId = which.max(cvEsts) #minimum location
    sdMax = sqrt(mean((R2sCv[maxId,]-cvEsts[maxId])^2)/(nrow(R2sCv)-1))
    cvId = cvEsts > (cvEsts[maxId] - sdMax)
    seId = which.max(cvId)
    cvR2 = cvEsts[seId]
    fullFits = lapply(lambdas, function(lam){ #Now the full fits with all lambdas
        gpls(data = data, glsSt = glsSt, xNames = xNames, outVar = outVar, corMat = corMat, lambda = lam)
    })
    coefs= sapply(fullFits, coef)
    outList = list("lambda" = lambdas, "cvm" = cvEsts, "cvsd" = sdMax, "coefs" = coef,
                   "lambda.min" = lambdas[maxId], "lambda.1se" = lamdas[seId])
}
