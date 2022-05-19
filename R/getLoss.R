#' Calculate the loss given predicted and observed values
#'
#' @param preds Matrix of predicted values
#' @param obs vector of observed values
#' @param loss a character vector indicating the loss type, see ?cv.pengls
#' @return the evaluated loss
getLoss = function(preds, obs, loss){
    SE = (preds-obs)^2
    switch(loss,
           "MSE" = SE,
           "R2" = 1-colMeans(SE)/var(obs)
    )
 }
