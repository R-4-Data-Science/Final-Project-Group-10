#' Identify Plausible Models Based on AIC and Stability
#'
#' `plausible_models()` filters a collection of candidate models (from a 
#' forward-search model path) by combining both AIC closeness and predictor 
#' stability scores. A model is considered "plausible" if:
#' 
#' * Its AIC is within `delta` of the minimum AIC among all explored models, and
#' * The average stability score of its predictors is at least `tau`.
#'
#' The function returns all models meeting both criteria.
#'
#' @param path_forest A list produced by `build_paths()` or similar functions, 
#'   containing at least:
#'   * `aic_by_model`: named list of AIC values (names are “model keys”), and  
#'   * `meta$predictors`: character vector of predictor names.
#'
#' @param path_stability A numeric vector produce by `stability()` or similar
#'   functions, giving stability scores for each predictor, in the same order as
#'   `path_forest$meta$predictors`. Stability values are typically in \[0, 1\], 
#'   where higher values indicate greater selection stability.
#'
#' @param delta Numeric. Maximum allowable AIC distance above the best model. 
#'   Models with `AIC <= min(AIC) + delta` are retained for further filtering. 
#'   Default is `2`.
#'
#' @param tau Numeric in \[0, 1\]. Minimum required average stability score for 
#'   a model to be considered plausible. Default is `0.6`.
#'
#' @details
#' The function operates in two stages:
#'
#' **Stage 1: AIC Filtering**  
#' All candidate models whose AIC lies within `delta` of the minimum observed AIC 
#' are collected.
#'
#' **Stage 2: Stability Filtering**  
#' For each model, the variables included are extracted. Their stability scores 
#' are averaged, and models with average stability less than `tau` are removed.
#'
#' The function returns only those models passing both filters.
#'
#' @return
#' A data frame with one row per plausible model, containing:
#'
#' * `model` – model key (string with variables joined by `"+"`)  
#' * `AIC` – numeric AIC value  
#' * `vars` – list column containing a character vector of predictors  
#' * `stability` – average stability score of the predictors in the model  
#'
#' @examples
#' \dontrun{
#'   # Suppose build_paths() was run:
#'   res <- build_paths(data = mtcars, response = "mpg",
#'                      predictors = colnames(mtcars)[-1])
#'
#'   # Suppose stability() was run:
#'   stab <- stability(data = mtcars, response = "mpg",
#'                      predictors = colnames(mtcars)[-1])
#'
#'   plausible_models(res, stab, delta = 2, tau = 0.6)
#' }
#'
#' @export
plausible_models <- function(path_forest, path_stability, delta=2, tau=0.6) {
  
  AICdf <- data.frame(
    model = names(path_forest$aic_by_model),
    AIC   = unlist(path_forest$aic_by_model),
    vars  = I(strsplit(names(path_forest$aic_by_model), "\\+"))
  )
  
  kept_models_index <- c()
  AICMin <- min(AICdf$AIC)
  
  for (i in 1:length(AICdf$model)) {
    if (AICdf$AIC[i] <= AICMin + delta) {
      kept_models_index <- c(kept_models_index, i)
    }
  }
  
  kept_models <- AICdf[kept_models_index, ]
  
  stability <- NA
  kept_models <- data.frame(kept_models, stability)
  
  remove.idx <- c()
  for (m in 1:length(kept_models$model)) {
    pidx <- which(path_forest$meta$predictors %in% kept_models$vars[[m]])
    
    kept_models$stability[m] <- mean(path_stability[pidx])
    if (kept_models$stability[m] < tau) {
      remove.idx <- c(remove.idx, m * -1)
    }
  }
  
  if (length(remove.idx) >= 1) {
    final_models <- kept_models[remove.idx, ]
  }
  else {
    final_models <- kept_models
  }
  
  return(final_models)
}
