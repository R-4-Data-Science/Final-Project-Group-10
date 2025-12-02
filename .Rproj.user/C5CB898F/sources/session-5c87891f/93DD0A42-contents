#' Resampling-based stability for multi-path selection
#'
#' For each resample of the data, this function runs [build_paths()] and
#' computes, for each predictor j, the proportion of models on that
#' resample that contain j. These per-resample proportions
#' z_j^(b) are then averaged over resamples to give stability scores
#' π_j = (1/B) Σ_b z_j^(b).
#'
#' @param x Predictor matrix or data.frame (n x p).
#' @param y Response vector of length n.
#' @param B Integer; number of resamples.
#' @param resample Character; type of resampling:
#'   either `"bootstrap"` or `"subsample"`.
#' @param m Integer; subsample size if `resample = "subsample"`.
#'   Ignored for bootstrap.
#' @param build_paths_fun Function used to run the multi-path search.
#'   Defaults to [build_paths()].
#' @param ... Additional arguments passed on to `build_paths_fun`
#'   (e.g., `family`, `K`, `eps`, `delta`, `L`).
#'
#' @return An object of class `"path_stability"` with components:
#' \describe{
#'   \item{pi}{Numeric vector of length p with stability scores π_j in [0,1].}
#'   \item{z}{B x p matrix with resample-wise proportions z_j^(b).}
#'   \item{meta}{List with resampling and build_paths metadata
#'     (B, resample, m, p, etc.).}
#' }
#'
#' @export
stability <- function(
    x, y,
    B = 50L,
    resample = c("bootstrap", "subsample"),
    m = NULL,
    build_paths_fun = build_paths,
    ...
) {
  resample <- match.arg(resample)
  
  x <- as.matrix(x)
  n <- nrow(x)
  p <- ncol(x)
  
  if (length(y) != n) {
    stop("Length of y (", length(y), ") must match nrow(x) (", n, ").")
  }
  
  # default subsample size if needed
  if (resample == "subsample" && is.null(m)) {
    m <- floor(sqrt(n))
  }
  
  colnames_x <- colnames(x)
  if (is.null(colnames_x)) {
    colnames_x <- paste0("x", seq_len(p))
  }
  
  z_mat <- matrix(NA_real_, nrow = B, ncol = p)
  colnames(z_mat) <- colnames_x
  
  for (b in seq_len(B)) {
    # draw indices for this resample
    idx <- if (resample == "bootstrap") {
      sample.int(n, n, replace = TRUE)
    } else {
      sample.int(n, m, replace = FALSE)
    }
    
    x_b <- x[idx, , drop = FALSE]
    y_b <- y[idx]
    
    # run multi-path search on this resample
    forest <- build_paths_fun(x_b, y_b, ...)
    
    # each model = vector of column indices
    models <- models_from_forest(forest)
    
    if (length(models) == 0L) {
      # no models: everything is zero for this resample
      z_mat[b, ] <- 0
    } else {
      M <- model_indicator_matrix(models, p = p)
      # proportion of models containing each predictor j on resample b
      z_mat[b, ] <- colMeans(M)
    }
  }
  
  pi <- colMeans(z_mat)
  
  out <- list(
    pi   = pi,
    z    = z_mat,
    meta = list(
      B        = B,
      resample = resample,
      m        = if (resample == "subsample") m else n,
      p        = p,
      colnames = colnames_x,
      build_paths_fun = deparse(substitute(build_paths_fun)),
      build_paths_args = list(...)
    )
  )
  
  class(out) <- c("path_stability", class(out))
  out
}

# Helper: extract list of models (each an integer vector of column indices)
# from a path_forest returned by build_paths().
models_from_forest <- function(forest) {
  frontiers <- forest$frontiers
  
  # Case 1: list of data.frames with a "vars" list-column
  if (is.list(frontiers) && all(vapply(frontiers, is.data.frame, logical(1)))) {
    models <- unlist(
      lapply(frontiers, function(df) df$vars),
      recursive = FALSE
    )
  } else if (is.list(frontiers) && all(vapply(frontiers, is.list, logical(1)))) {
    # Case 2: list of lists-of-models
    models <- unlist(frontiers, recursive = FALSE)
  } else {
    stop("models_from_forest(): don't know how to extract models from 'frontiers'. ",
         "Adapt this helper to your build_paths() output.")
  }
  
  # drop duplicates: sort each model, then unique on string key
  models <- lapply(models, sort)
  keys   <- vapply(models, paste, collapse = ",", FUN.VALUE = character(1))
  models[!duplicated(keys)]
}

# Helper: given a list of models, build an indicator matrix:
# rows = models, cols = predictors 1..p
model_indicator_matrix <- function(models, p) {
  n_models <- length(models)
  M <- matrix(0L, nrow = n_models, ncol = p)
  
  if (n_models == 0L) return(M)
  
  for (k in seq_len(n_models)) {
    idx <- models[[k]]
    if (length(idx) > 0L) {
      M[k, idx] <- 1L
    }
  }
  
  M
}

