# Final-Project-Group-10

This repository contains the source code for the `Stepwise` R package.  
The package implements a stepwise model selection procedure based on
multi-path AIC searches with a resampling stability assessment and
the identification of plausible models for linear and logistic regression.

The typical workflow is:

1. Build a multi-path forest of models using `forward_path_selection()`.
2. Estimate selection stability via resampling using `stability()`.
3. Select a set of plausible models using `plausible_models()`.

## Installation Procedure

```r
install.packages("remotes")
remotes::install_github("R-4-Data-Science/Final-Project-Group-10/Stepwise")
library(Stepwise)
```
## Example with linear regression:
```r
set.seed(1)
n <- 120
p <- 8

X <- matrix(rnorm(n * p), n, p)
beta <- c(2, -1.5, 0, 0, 1, rep(0, p - 5))
y <- as.numeric(X %*% beta + rnorm(n, sd = 1))

colnames(X) <- paste0("x", 1:p)
df <- data.frame(y = y, X)
predictors <- setNames(paste0("x", 1:p), paste0("x", 1:p))

# 1) Multi-path forward selection
forest <- build_paths(
  data       = df,
  response   = "y",
  predictors = names(predictors),
  # family    = "gaussian",   
  K          = min(p, 10),
  eps        = 1e-6,
  delta      = 1,
  L          = 50,
  verbose    = FALSE
)

# 2) Resampling-based stability
pi <- stability(
  data       = df,
  response   = "y",
  predictors = names(predictors),
  B          = 50,
  resample   = "bootstrap",
  m          = NULL,       
  delta      = 1,
  eps        = 1e-6,
  L          = 50,
  K          = min(p, 10),
  verbose    = FALSE
)

# 3) Plausible models filtered by AIC window + average stability
plaus_lin <- plausible_models(
  forest = forest,
  pi     = pi,
  Delta  = 2,   
  tau    = 0.6   
)

plaus_lin
```
##Example with logistic regression:
```r
set.seed(2)
n <- 200
p <- 6

Xb <- matrix(rnorm(n * p), n, p)
linpred <- 1.2 * Xb[, 1] - 1.0 * Xb[, 2] + 0.8 * Xb[, 5]
prob    <- 1 / (1 + exp(-linpred))
ybin    <- rbinom(n, size = 1, prob = prob)

colnames(Xb) <- paste0("x", 1:p)
dfb <- data.frame(y = ybin, Xb)
predictors_bin <- paste0("x", 1:p)

# 1) Multi-path forward selection (binomial family)
forest_bin <- build_paths(
  data       = dfb,
  response   = "y",
  predictors = predictors_bin,
  # family    = "binomial",  
  K          = min(p, 10),
  eps        = 1e-6,
  delta      = 1,
  L          = 50,
  verbose    = FALSE
)

# 2) Resampling-based stability for logistic model
pi_bin <- stability(
  data       = dfb,
  response   = "y",
  predictors = predictors_bin,
  B          = 50,
  resample   = "bootstrap",
  m          = NULL,
  delta      = 1,
  eps        = 1e-6,
  L          = 50,
  K          = min(p, 10),
  verbose    = FALSE
)

# 3) Plausible models for the logistic case
plaus_logit <- plausible_models(
  forest = forest_bin,
  pi     = pi_bin,
  Delta  = 2,
  tau    = 0.6
)

plaus_logit
```
