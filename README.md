
This repository demonstrates the use of the *pengls* package for
high-dimensional data with spatial or temporal autocorrelation. It
consists of an iterative loop around the *nlme* and *glmnet* packages.
Currently, only continuous outcomes and $R^2$ as performance measure are
implemented.

# Installation instuctions

The *pengls* package is available from BioConductor, and can be
installed as follows:

``` r
library(BiocManager)
install("pengls")
```

Once installed, it can be loaded and version info printed.

``` r
suppressPackageStartupMessages(library(pengls))
cat("pengls package version", as.character(packageVersion("pengls")), "\n")
```

    ## pengls package version 1.5.0

# Illustration

We first create a toy dataset with spatial coordinates.

``` r
library(nlme)
n <- 75 #Sample size
p <- 100 #Number of features
g <- 10 #Size of the grid
#Generate grid
Grid <- expand.grid("x" = seq_len(g), "y" = seq_len(g))
# Sample points from grid without replacement
GridSample <- Grid[sample(nrow(Grid), n, replace = FALSE),]
#Generate outcome and regressors
b <- matrix(rnorm(p*n), n , p)
a <- rnorm(n, mean = b %*% rbinom(p, size = 1, p = 0.25), sd = 0.1) #25% signal
#Compile to a matrix
df <- data.frame("a" = a, "b" = b, GridSample)
```

The *pengls* method requires prespecification of a functional form for
the autocorrelation. This is done through the *corStruct* objects
defined by the *nlme* package. We specify a correlation decaying as a
Gaussian curve with distance, and with a nugget parameter. The nugget
parameter is a proportion that indicates how much of the correlation
structure explained by independent errors; the rest is attributed to
spatial autocorrelation. The starting values are chosen as reasonable
guesses; they will be overwritten in the fitting process.

``` r
# Define the correlation structure (see ?nlme::gls), with initial nugget 0.5 and range 5
corStruct <- corGaus(form = ~ x + y, nugget = TRUE, value = c("range" = 5, "nugget" = 0.5))
```

Finally the model is fitted with a single outcome variable and large
number of regressors, with the chosen covariance structure and for a
prespecified penalty parameter $\lambda=0.2$.

``` r
#Fit the pengls model, for simplicity for a simple lambda
penglsFit <- pengls(data = df, outVar = "a", xNames = grep(names(df), pattern = "b", value =TRUE),
glsSt <- corStruct, lambda = 0.2, verbose = TRUE)
```

    ## Starting iterations...
    ## Iteration 1 
    ## Iteration 2 
    ## Iteration 3

Standard extraction functions like print(), coef() and predict() are
defined for the new “pengls” object.

``` r
penglsFit
```

    ## pengls model with correlation structure: corGaus 
    ##  and 47 non-zero coefficients

``` r
penglsCoef <- coef(penglsFit)
penglsPred <- predict(penglsFit)
```
