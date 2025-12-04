# 3.3
# https://chatgpt.com/share/6931ad54-5900-8011-a70d-63be4254d81c


plausible_models <- function(path_forest, path_stability, delta=2, tau=0.6) {
  
  AICdf <- data.frame(model=names(path_forest$aic_by_model), AIC=unlist(path_forest$aic_by_model), vars=I(strsplit(names(path_forest$aic_by_model), "\\+")))
  kept_models_index <- c()
  AICMin <- min(AICdf$AIC)
  
  for (i in 1:length(AICdf$model)) {
    if (AICdf$AIC[i] <= AICMin + delta) {
      kept_models_index <- c(kept_models_index, i)
    }
  }
  
  kept_models <- AICdf[kept_models_index,]
  
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
  
  final_models <- kept_models[remove.idx,]
  
  return(final_models)
}
