#' Stability estimation with resampling
#'
#' Estimates selection stability by running \code{build_paths()} on resampled
#' data and recording how often each predictor appears in the models.
#'
#' @param data Data frame containing the response and predictors.
#' @param response Name of the response variable (character).
#' @param predictors Character vector of predictor column names in \code{data}.
#' @param B Number of resamples (bootstrap samples or subsamples).
#' @param resample Resampling scheme: \code{"bootstrap"} or \code{"subsample"}.
#' @param m Subsample size when \code{resample = "subsample"}; defaults to
#'   \code{floor(0.8 * nrow(data))} if \code{NULL}.
#' @param delta Tuning parameter passed to \code{build_paths()} controlling the
#'   allowable AIC distance from the best child.
#' @param eps Tuning parameter passed to \code{build_paths()} for the minimum
#'   AIC improvement over the parent.
#' @param L Maximum number of models kept at each step in \code{build_paths()}.
#' @param K Maximum number of forward steps (defaults to \code{length(predictors)}).
#' @param verbose Logical; if \code{TRUE}, prints progress over resamples.
#'
#' @return
#' A named numeric vector of length \code{length(predictors)}, with values in
#' \code{[0, 1]} giving the estimated stability of each predictor.
#'
#'
#' @importFrom stats setNames
#' @examples
#' \dontrun{
#' set.seed(1)
#' dat <- data.frame(
#'   y  = rnorm(50),
#'   x1 = rnorm(50),
#'   x2 = rnorm(50),
#'   x3 = rnorm(50)
#' )
#'
#' stability(
#'   data       = dat,
#'   response   = "y",
#'   predictors = c("x1", "x2", "x3"),
#'   B          = 10,
#'   resample   = "bootstrap"
#' )
#' }
#' @export
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

  ## stability vector π_j = (1/B) * sum_b z_j^(b)
  pi_sum / B
}

