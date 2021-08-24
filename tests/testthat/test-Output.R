context("Output components")
library(nlme)
n = 50 #Sample size
p = 100 #Number of features
g = 10 #Size of the grid
#Generate grid
Grid = expand.grid("x" = seq_len(g), "y" = seq_len(g))
# Sample points from grid without replacement
GridSample = Grid[sample(nrow(Grid), n, replace = FALSE),]
#Generate outcome and regressors
a = rnorm(n)
b = matrix(rnorm(p*n), n , p)
#Compile to a matrix
df = data.frame("a" = a, "b" = b, GridSample)
# Define the correlation structure (see ?nlme::gls), with initial nugget 0.5 and range 5
corStruct = corGaus(form = ~ x + y, nugget = TRUE, value = c("range" = 5, "nugget" = 0.5))
#Fit the gpls model, for simplicity for a simple lambda
gplsFit = gpls(data = df, outVar = "a", xNames = grep(names(df), pattern = "b", value =TRUE),
glsSt = corStruct, nfolds = 5)
test_that("Function returns list with correct elements", {
    expect_s3_class(gplsFit, "gpls")
    expect_output(print(gplsFit))
    expect_s4_class(coef(gplsFit), "dgCMatrix")
}
)
