#' Peform cross-validation pengls
#'
#' @inheritParams pengls
#' @param lambdas an optional lambda sequence
#' @param transFun a transformation function to apply to predictions and outcome in the cross-validation
#' @param transFunArgs Additional arguments passed onto transFun
#'
#' @return A list with components
#' \item{lambda}{The series of lambdas}
#' \item{cvm}{The vector of mean R2's}
#' \item{cvsd}{The standard error of R2 at the maximum}
#' \item{cvOpt}{The R2 according to the 1 standard error rule}
#' \item{coefs}{The matrix of coefficients for every lambda value}
#' \item{lambda.min}{Lambda value with maximal R2}
#' \item{lambda.1se}{Smallest lambda value within 1 standard error from the maximum}
#' \item{foldid}{The folds}
#' \item{glsSt}{The nlme correlation object}
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
#' corStruct = corGaus(form = ~ x + y, nugget = TRUE, value = c("range" = 5, "nugget" = 0.5))
#' #Fit the pengls model, for simplicity for a simple lambda
#' register(MulticoreParam(3)) #Prepare multithereading
#' penglsFitCV = cv.pengls(data = df, outVar = "a",
#' xNames = grep(names(df), pattern = "b", value =TRUE),
#' glsSt = corStruct, nfolds = 5)
#' penglsFitCV$lambda.1se #Lambda for 1 standard error rule
#' penglsFitCV$cvOpt #Corresponding R2
#' coef(penglsFitCV)
#' penglsFitCV$foldid #The folds used
cv.pengls = function(data, glsSt, xNames, outVar, corMat, nfolds, foldid,
                   cvType = "blocked", lambdas, transFun = "identity",
                   transFunArgs = list(), ...){
    if(missing(foldid))
        foldid <- makeFolds(nfolds, data, cvType)
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
    R2sCv <- vapply(FUN.VALUE = lambdas, seq_along(fits), function(i){
        id = foldid==uniqueFolds[i]
        coefs = vapply(FUN.VALUE = numeric(length(xNames)+1), fits[[i]], function(x) as.vector(coef(x)))
        preds = apply(cbind(1, as.matrix(data[id,xNames, drop = FALSE])) %*% coefs, 2,
                      function(x) do.call(what = transFun, args = c(transFunArgs, list("x" =x))))
        testDat = do.call(what = transFun, args = c(transFunArgs, list("x"= data[id,outVar])))
        1-colMeans((preds-testDat)^2)/mean((testDat-mean(testDat))^2)
    })
    #Use lambda within 1se from minimum error
    cvEsts <- rowMeans(R2sCv)
    maxId <- which.max(cvEsts) #minimum location
    sdMax <- sqrt(mean((R2sCv[maxId,]-cvEsts[maxId])^2)/(nrow(R2sCv)-1))
    cvId <- cvEsts > (cvEsts[maxId] - sdMax)
    seId <- which.max(cvId)

    cvR2 <- cvEsts[seId]
    fullFits <- bplapply(lambdas, function(lam){ #Now the full fits with all lambdas
        pengls(data = data, glsSt = glsSt, xNames = xNames, outVar = outVar, lambda = lam)
    })
    coefs <- vapply(FUN.VALUE = numeric(length(xNames)+1), fullFits, function(x) as.vector(coef(x)))
    outList <- list("lambda" = lambdas, "cvm" = cvEsts, "cvsd" = sdMax,
                   "cvOpt" = cvEsts[seId], "coefs" = coefs,
                   "lambda.min" = lambdas[maxId], "lambda.1se" = lambdas[seId],
                   "foldid" = foldid, "glsSt" = glsSt)
    class(outList) <- "cv.pengls"
    return(outList)
}
