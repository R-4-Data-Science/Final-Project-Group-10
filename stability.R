## 3.2  Stability estimation with resampling
stability <- function(
    data,
    response,
    predictors,
    B        = 50,
    resample = c("bootstrap", "subsample"),
    m        = NULL,                  
    delta    = 2,
    eps      = 0.5,
    L        = 20,
    K        = length(predictors),
    verbose  = TRUE
) {
  resample <- match.arg(resample)
  n        <- nrow(data)
  p        <- length(predictors)
  
  # accumulator for sum_b z_j^(b)
  pi_sum <- setNames(numeric(p), predictors)
  
  for (b in seq_len(B)) {
    
    ## 1. draw resample indices
    idx <- if (resample == "bootstrap") {
      sample.int(n, n, replace = TRUE)
    } else {
      if (is.null(m)) m <- floor(0.8 * n)
      sample.int(n, m, replace = FALSE)
    }
    
    data_b <- data[idx, , drop = FALSE]
    
    ## 2. run multi-path search 
    paths_b <- build_paths(
      data       = data_b,
      response   = response,
      predictors = predictors,
      delta      = delta,
      eps        = eps,
      L          = L,
      K          = K,
      verbose    = FALSE
    )$path_forest
    
    ## 3. collect all models across steps, drop empty model
    models_b <- unlist(lapply(paths_b, function(df) df$model), recursive = FALSE)
    models_b <- Filter(function(m) length(m) > 0, models_b)
    
    ## 4. z_j^(b) = proportion of models containing predictor j
    z_b <- numeric(p)
    if (length(models_b) > 0) {
      for (j in seq_along(predictors)) {
        v <- predictors[j]
        z_b[j] <- mean(vapply(models_b, function(m) v %in% m, logical(1)))
      }
    }
    
    ## accumulate
    pi_sum <- pi_sum + z_b
    if (verbose) message("Resample ", b, " of ", B)
  }
  
  ## 5. stability vector π_j = (1/B) * sum_b z_j^(b)
  pi_sum / B
}

## quick test w mtcars

data(mtcars)
response   <- "mpg"
predictors <- setdiff(colnames(mtcars), response)

set.seed(123)
pi_mtcars <- stability(
  data       = mtcars,
  response   = response,
  predictors = predictors,
  B          = 20,
  resample   = "bootstrap",
  delta      = 2,
  eps        = 0.5,
  L          = 20,
  verbose    = FALSE
)

pi_mtcars
