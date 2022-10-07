context("Output components")
library(nlme)
n <- 50 #Sample size
p <- 100 #Number of features
g <- 10 #Size of the grid
#Generate grid
Grid <- expand.grid("x" = seq_len(g), "y" = seq_len(g))
# Sample points from grid without replacement
GridSample <- Grid[sample(nrow(Grid), n, replace = FALSE),]
#Generate outcome and regressors
a <- rnorm(n)
b <- matrix(rnorm(p*n), n , p)
#Compile to a matrix
df <- data.frame("out" = a, "b" = b, GridSample)
# Define the correlation structure (see ?nlme::gls), with initial nugget 0.5 and range 5
corStruct <- corGaus(form = ~ x + y, nugget = TRUE, value = initValues <- c("range" = 5, "nugget" = 0.5))
#Fit the pengls model, for simplicity for a simple lambda
penglsFit <- pengls(data = df, outVar = "out", xNames = grep(names(df), pattern = "b", value =TRUE),
glsSt <- corStruct, nfolds = 5)
test_that("pengls function returns list with correct elements", {
    expect_s3_class(penglsFit, "pengls")
    expect_type(predict(penglsFit), "double")
    expect_output(print(penglsFit))
    expect_type(coef(penglsFit), "double")
}
)
penglsFitCv <- cv.pengls(data = df, outVar = "out", xNames = grep(names(df), pattern = "b", value =TRUE),
               glsSt = corStruct, nfolds = 5)
test_that("cv.pengls function returns list with correct elements", {
    expect_s3_class(penglsFitCv, "cv.pengls")
    expect_output(print(penglsFitCv))
    expect_type(coef(penglsFitCv), "double")
    expect_type(predict(penglsFitCv), "double")
}
)
test_that("Parameters of correlation matrix are being estimated", {
    expect_false(any(penglsFit$gls$Coef==initValues))
}
)
