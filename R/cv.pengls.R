#' Peform cross-validation pengls
#'
#' @inheritParams pengls
#' @param lambdas an optional lambda sequence
#' @param transFun a transformation function to apply to predictions and outcome in the cross-validation
#' @param transFunArgs Additional arguments passed onto transFun
#' @param cvType A character vector defining the type of cross-validation.
#' Either "random" or "blocked", ignored if foldid is provided
#' @param loss a character vector, currently either 'R2' or 'MSE' indicating the
#'  loss function (although R2 is not a proper loss...)
#'
#' @return A list with components
#' \item{lambda}{The series of lambdas}
#' \item{cvm}{The vector of mean R2's}
#' \item{cvsd}{The standard error of R2 at the maximum}
#' \item{cvOpt}{The R2 according to the 1 standard error rule}
#' \item{coefs}{The matrix of coefficients for every lambda value}
#' \item{bestFit}{The best fitting pengls model according to the 1 standard error rule}
#' \item{lambda.min}{Lambda value with maximal R2}
#' \item{lambda.1se}{Smallest lambda value within 1 standard error from the maximum}
#' \item{foldid}{The folds}
#' \item{glsSt}{The nlme correlation object}
#' \item{loss}{The loss function used}
#' @export
#' @importFrom BiocParallel bplapply
#'
#' @examples
#' library(nlme)
#' library(BiocParallel)
#' n <- 20 #Sample size
#' p <- 50 #Number of features
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
#' corStruct = corGaus(form = ~ x + y, nugget = TRUE,
#' value = c("range" = 5, "nugget" = 0.5))
#' #Fit the pengls model, for simplicity for a simple lambda
#' register(MulticoreParam(3)) #Prepare multithereading
#' penglsFitCV = cv.pengls(data = df, outVar = "a", xNames = grep(names(df),
#' pattern = "b", value = TRUE),
#' glsSt = corStruct, nfolds = 5)
#' penglsFitCV$lambda.1se #Lambda for 1 standard error rule
#' penglsFitCV$cvOpt #Corresponding R2
#' coef(penglsFitCV)
#' penglsFitCV$foldid #The folds used
#' #With MSE as loss function
#' penglsFitCVmse = cv.pengls(data = df, outVar = "a",
#' xNames = grep(names(df), pattern = "b", value =TRUE),
#' glsSt = corStruct, nfolds = 5, loss = "MSE")
#' penglsFitCVmse$lambda.1se #Lambda for 1 standard error rule
#' penglsFitCVmse$cvOpt #Corresponding MSE
#' coef(penglsFitCVmse)
#' predict(penglsFitCVmse)
cv.pengls = function(data, glsSt, xNames, outVar, corMat, nfolds, foldid, scale = TRUE, center = TRUE,
                   cvType = "blocked", lambdas, transFun = "identity",
                   transFunArgs = list(), loss = c("R2", "MSE"), ...){
    loss = match.arg(loss)
    if(missing(foldid))
        foldid <- makeFolds(nfolds, data, cvType)
    data[, xNames] = scale(data[, xNames], scale = scale, center = center) #Center and scale
    if(missing(lambdas)){
        naieveFit <- cv.glmnet(x = as.matrix(data[,xNames]),
                              y = data[,outVar], foldid = foldid, ...)
        #The naive fit provides a lambda sequence
        lambdas <- naieveFit$lambda
    }
    uniqueFolds <- unique(foldid)
    fits <- bplapply(uniqueFolds, function(i){
        id = foldid!=i #Leave out test fold
        lapply(lambdas, function(lam){
                pengls(data = data[id,], glsSt = glsSt, xNames = xNames,
                     outVar = outVar, lambda = lam)
        })
    })
    #Find R2 on left out folds
    LossCv0 <- lapply(seq_along(fits), function(i){
        id = foldid==uniqueFolds[i]
        coefs = vapply(FUN.VALUE = numeric(length(xNames)+1), fits[[i]], function(x) as.vector(coef(x)))
        preds = apply(cbind(1, as.matrix(data[id,xNames, drop = FALSE])) %*% coefs, 2,
                      function(x) do.call(what = transFun, args = c(transFunArgs, list("x" =x))))
        testDat = do.call(what = transFun, args = c(transFunArgs, list("x"= data[id,outVar])))
        getLoss(preds, testDat, loss)
    })
    LossCv = switch(loss,
                    "R2" = Reduce(LossCv0, f = cbind),
                    "MSE" = t(Reduce(LossCv0, f = rbind)))
    #Use lambda within 1se from minimum error
    cvEsts <- rowMeans(LossCv)
    whichFun = switch(loss, "MSE" = which.min, "R2" = which.max)
    maxId <- whichFun(cvEsts) #minimum location
    sdMax <- sqrt(mean((LossCv[maxId,]-cvEsts[maxId])^2)/(nrow(LossCv)-1))
    cvId <- switch(loss, "MSE" = cvEsts < (cvEsts[maxId] + sdMax),
                   "R2" = cvEsts > (cvEsts[maxId] - sdMax))
    seId <- which.max(cvId)

    cvR2 <- cvEsts[seId]
    fullFits <- bplapply(lambdas, function(lam){ #Now the full fits with all lambdas
        pengls(data = data, glsSt = glsSt, xNames = xNames, outVar = outVar, lambda = lam)
    })
    coefs <- vapply(FUN.VALUE = numeric(length(xNames)+1), fullFits, function(x) as.vector(coef(x)))
    outList <- list("lambda" = lambdas, "cvm" = cvEsts, "cvsd" = sdMax,
                   "cvOpt" = cvEsts[seId], "coefs" = coefs, "bestFit" = fullFits[[seId]],
                   "lambda.min" = lambdas[maxId], "lambda.1se" = lambdas[seId],
                   "foldid" = foldid, "glsSt" = glsSt, "loss" = loss)
    class(outList) <- "cv.pengls"
    return(outList)
}
