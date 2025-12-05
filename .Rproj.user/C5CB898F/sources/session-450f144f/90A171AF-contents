#' Build forward selection paths
#'
#' Performs forward-path selection starting from the empty model, adding one variable at a time,
#' keeping children within delta of the best AIC and exceeding minimum improvement eps.
#' If no children meet eps, the best child is always kept to allow multistep selection.
#' Deduplicates and prunes to a maximum L models per step.
#'
#' @param data A data frame containing the response and predictors
#' @param response Name of the response variable (string)
#' @param predictors Vector of predictor variable names (strings)
#' @param delta Numeric. Allowable AIC distance from the parent’s best child (default 2)
#' @param eps Numeric. Minimum required AIC improvement over parent (default 0.5)
#' @param L Integer. Maximum number of models kept per step (default 20)
#' @param K Integer. Maximum number of steps / model size (default = length(predictors))
#' @param verbose Logical. Whether to print progress (default TRUE)
#' @param debug Logical. Whether to print per-step debug info (default FALSE)
#'
#' @return A list with:
#'   - path_forest: list of data frames with models and AIC per step
#'   - aic_by_model: named vector of AIC for each model
#'   - meta: list with parameters and info
#'
#' @examples
#' set.seed(1)
#' n <- 120; p <- 8
#' X <- matrix(rnorm(n*p), n, p)
#' colnames(X) <- paste0("x", 1:p)
#' beta <- c(2, -1.5, 0, 0, 1, rep(0, p-5))
#' y <- as.numeric(X %*% beta + rnorm(n))  # numeric vector
#' dat <- data.frame(y=y, X)
#' preds <- colnames(dat)[-1]
#'
#' # Run forward selection with multiple steps
#' res <- build_paths(dat, "y", preds, delta=2, eps=0.5, L=5, K=4, verbose=TRUE, debug=FALSE)
#'
#' # Inspect frontiers at each step
#' for(i in seq_along(res$path_forest)) {
#'   cat("Step", i-1, "frontier models:\n")
#'   print(res$path_forest[[i]])
#'   cat("\n")
#' }
#'
#' # Access metadata
#' res$meta
build_paths <- function(
    data,
    response,
    predictors,
    delta = 2,
    eps = 0.5,
    L = 20,
    K = length(predictors),
    verbose = TRUE,
    debug = FALSE
) {
  stopifnot(response %in% colnames(data))

  model_key <- function(vars) paste(sort(vars), collapse = "+")
  model_aic <- function(vars) {
    f <- as.formula(paste(response, "~", ifelse(length(vars)==0, "1", paste(vars, collapse="+"))))
    AIC(lm(f, data = data)) #for logical change lm() → glm(..., family = binomial) and ensure response is a factor 
  }

  aic_by_model <- list()
  path_forest <- list()

  # Step 0: empty model
  empty <- character(0)
  null_aic <- model_aic(empty)
  aic_by_model[[model_key(empty)]] <- null_aic
  path_forest[[1]] <- data.frame(
    model = I(list(empty)),
    AIC = null_aic
  )

  step <- 1
  repeat {
    if (step > K) break  # stop if step limit reached

    parents <- lapply(path_forest[[step]]$model, unlist)
    if (verbose) message("Step ", step, ": ", length(parents), " parent models")

    children_all <- list()
    child_keys <- character()

    for (parent in parents) {
      parent_k <- model_key(parent)
      parent_aic <- aic_by_model[[parent_k]]

      remaining <- setdiff(predictors, parent)
      if (length(remaining) == 0) next

      if (debug) {
        cat("Parent model:", ifelse(length(parent)==0, "(empty)", paste(parent, collapse="+")),
            "AIC =", parent_aic, "\n")
      }

      child_info <- lapply(remaining, function(v) {
        vars <- c(parent, v)
        key <- model_key(vars)
        if (is.null(aic_by_model[[key]])) {
          aic_by_model[[key]] <<- model_aic(vars)
        }
        list(vars = vars, key = key, AIC = aic_by_model[[key]])
      })

      child_AICs <- sapply(child_info, `[[`, "AIC")
      best_child <- min(child_AICs)

      # Keep children within delta and improvement >= eps
      keep_idx <- which((child_AICs - best_child <= delta) &
                          (parent_aic - child_AICs >= eps))

      # FORCE keep the best child if none satisfy eps
      if (length(keep_idx) == 0) {
        keep_idx <- which.min(child_AICs)
      }

      if (debug) {
        for (i in seq_along(child_info)) {
          vars <- child_info[[i]]$vars
          aic <- child_info[[i]]$AIC
          improvement <- parent_aic - aic
          keep <- ifelse(i %in% keep_idx, "KEEP", "REJECT")
          cat("  Child:", paste(vars, collapse="+"), "AIC =", aic,
              "improvement =", improvement, keep, "\n")
        }
      }

      kept <- child_info[keep_idx]
      children_all <- c(children_all, lapply(kept, `[[`, "vars"))
      child_keys <- c(child_keys, sapply(kept, `[[`, "key"))
    }

    if (length(child_keys) == 0) {
      if (verbose) message("No children left. Stopping.")
      break
    }

    # Deduplicate
    unique_idx <- !duplicated(child_keys)
    frontier <- children_all[unique_idx]

    # Prune if > L
    if (length(frontier) > L) {
      AICs <- sapply(frontier, function(v) aic_by_model[[model_key(v)]])
      frontier <- frontier[order(AICs)[1:L]]
    }

    # Store frontier
    df_frontier <- data.frame(
      model = I(frontier),
      AIC = sapply(frontier, function(v) aic_by_model[[model_key(v)]])
    )
    path_forest[[step + 1]] <- df_frontier

    step <- step + 1
  }

  list(
    path_forest = path_forest,
    aic_by_model = aic_by_model,
    meta = list(
      response = response,
      predictors = predictors,
      delta = delta,
      eps = eps,
      L = L,
      K = K
    )
  )
}
