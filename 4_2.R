## 4.2  Stability Algorithm

stability_glm <- function(
    X,
    y,
    family   = c("gaussian", "binomial"),
    B        = 50,
    resample = c("bootstrap", "subsample"),
    m        = NULL,      
    delta    = 2,
    eps      = 0.5,
    L        = 20,
    K        = 5,         
    verbose  = TRUE
) {
  family   <- match.arg(family)
  resample <- match.arg(resample)
  
  X <- as.data.frame(X)
  n <- nrow(X)
  
  predictors <- colnames(X)
  p          <- length(predictors)
  
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
    
    X_b <- X[idx, , drop = FALSE]
    y_b <- y[idx]
    
    tree_b <- algorithm_forward_tree(
      X       = X_b,
      y       = y_b,
      family  = family,
      K       = K,
      eps     = eps,
      delta   = delta,
      L       = L,
      verbose = FALSE
    )
    
    models_by_step <- tree_b$models_by_step  # list of lists of character vectors
    
    ##collect all models across steps, drop empty model
    models_b <- unlist(models_by_step, recursive = FALSE)
    models_b <- Filter(function(m) length(m) > 0, models_b)
    
    ##z_j^(b) = proportion of models containing predictor j
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
  
  ##stability vector π_j = (1/B) * sum_b z_j^(b)
  pi_sum / B
}
