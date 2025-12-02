#3.1

 
forward_path_selection <- function(
    data,
    response,
    predictors,
    delta = 2,    # allowable AIC distance from parent’s best child
    eps = 0.5,    # required AIC improvement over parent
    L = 20,       # maximum number of models kept per step
    verbose = TRUE
) {
  
  # ---- helpers -------------------------------------------------------------
  
  model_key <- function(vars) paste(sort(vars), collapse = "+")  # unique model key
  
  model_aic <- function(vars) {                                 # fit model + compute AIC
    f <- as.formula(
      paste(response, "~",
            ifelse(length(vars) == 0, "1", paste(vars, collapse = "+")))
    )
    AIC(lm(f, data = data))
  }
  
  # ---- storage -------------------------------------------------------------
  
  aic_by_model <- list()
  path_forest  <- list()
  
  # Step 0: empty model
  empty <- character(0)
  path_forest[[1]] <- list(empty)
  aic_by_model[[model_key(empty)]] <- model_aic(empty)
  
  # ---- forward selection ---------------------------------------------------
  
  step <- 1
  repeat {
    
    parents <- path_forest[[step]]
    if (verbose) message("Step ", step, ": ", length(parents), " parents.")
    
    children_all <- list()
    child_keys <- character()
    
    for (parent in parents) {
      
      parent_k   <- model_key(parent)
      parent_aic <- aic_by_model[[parent_k]]
      
      remaining <- setdiff(predictors, parent)
      if (length(remaining) == 0) next
      
      # Generate all child models (add exactly one variable)
      child_info <- lapply(remaining, function(v) {
        vars <- c(parent, v)
        key  <- model_key(vars)
        
        if (is.null(aic_by_model[[key]])) {
          aic_by_model[[key]] <<- model_aic(vars)
        }
        
        list(vars = vars, key = key, AIC = aic_by_model[[key]])
      })
      
      child_AICs <- sapply(child_info, `[[`, "AIC")
      best_child <- min(child_AICs)
      
      # keep children improving AIC enough and within delta
      keep_idx <- which((child_AICs - best_child <= delta) &
                          (parent_aic - child_AICs >= eps))
      
      if (length(keep_idx) > 0) {
        kept <- child_info[keep_idx]
        children_all <- c(children_all, lapply(kept, `[[`, "vars"))
        child_keys   <- c(child_keys, sapply(kept, `[[`, "key"))
      }
    }
    
    if (length(child_keys) == 0) {
      if (verbose) message("No AIC-improving children. Stopping.")
      break
    }
    
    # Remove duplicates
    unique_idx <- !duplicated(child_keys)
    frontier <- children_all[unique_idx]
    
    # prune if needed
    if (length(frontier) > L) {
      if (verbose) message("Pruning from ", length(frontier), " → best ", L, " by AIC.")
      AICs <- sapply(frontier, function(v) aic_by_model[[model_key(v)]])
      frontier <- frontier[order(AICs)[1:L]]
    }
    
    path_forest[[step + 1]] <- frontier
    step <- step + 1
  }
  
  # ---- output --------------------------------------------------------------
  
  list(
    path_forest  = path_forest,
    aic_by_model = aic_by_model,
    meta = list(
      response   = response,
      predictors = predictors,
      delta = delta,
      eps   = eps,
      L     = L
    )
  )
}
