#3.1

build_paths <- function(y, X) {
  p <- ncol(X)
  colnames(X) <- paste0("x", 1:p)
  
  # Put everything in a data frame for lm()
  data <- data.frame(y = y, X)
  
  null_model <- lm(y ~ 1, data = data)
  null_aic <- AIC(null_model)
  
  parents <- list(c())  # start with empty model
  step <- 1
  
  repeat {
    cat("Step", step, ":", length(parents), "parents.\n")
    children <- list()
    child_aic <- c()
    
    # Loop over parents
    for (parent in parents) {
      
      parent_formula <- as.formula(
        paste("y ~", ifelse(length(parent)==0, "1", paste(parent, collapse="+")))
      )
      parent_model <- lm(parent_formula, data = data)
      parent_aic <- AIC(parent_model)
      
      unused <- setdiff(colnames(X), parent)
      
      # Try adding each unused variable
      for (u in unused) {
        new_terms <- c(parent, u)
        fmla <- as.formula(paste("y ~", paste(new_terms, collapse="+")))
        mod <- lm(fmla, data = data)
        a <- AIC(mod)
        improvement <- parent_aic - a
        
        cat("DEBUG:", deparse(fmla), 
            " AIC=", round(a,4),
            " improvement=", round(improvement,4), "\n")
        
        children[[length(children)+1]] <- new_terms
        child_aic <- c(child_aic, a)
      }
    }
    
    # Evaluate whether children improve over parents
    parent_aics <- sapply(parents, function(pr) {
      f <- as.formula(
        paste("y ~", ifelse(length(pr)==0,"1",paste(pr,collapse="+")))
      )
      AIC(lm(f, data = data))
    })
    
    best_parent_aic <- min(parent_aics)
    
    improving <- which(child_aic < best_parent_aic)
    
    if (length(improving) == 0) {
      cat("No AIC-improving children. Stopping.\n")
      break
    }
    
    parents <- children[improving]
    step <- step + 1
  }
  
  return(parents)
}

#4.1

algorithm_forward_tree <- function(
    X, 
    y, 
    family = c("gaussian", "binomial"),
    K = 5,
    eps = 0.5,
    delta = 2,
    L = 20,
    verbose = TRUE
) {
  
  family <- match.arg(family)
  predictors <- colnames(X)
  
  # Helper: convert to unique model key
  model_key <- function(vars) paste(sort(vars), collapse = "+")
  
  # Helper: fit GLM and compute AIC
  model_aic <- function(vars) {
    df <- data.frame(y = y, X)
    if (length(vars) == 0) {
      f <- as.formula("y ~ 1")
    } else {
      f <- as.formula(paste("y ~", paste(vars, collapse = "+")))
    }
    AIC(glm(f, data = df, family = family))
  }
  
  # Storage
  aic_by_model <- list()
  model_tree   <- list()
  
  # Step 0 — empty model
  empty <- character(0)
  model_tree[[1]] <- list(empty)
  aic_by_model[[model_key(empty)]] <- model_aic(empty)
  
  # --------------------------------------------------------------------------
  # Main Loop
  # --------------------------------------------------------------------------
  
  for (k in 1:K) {
    
    parents <- model_tree[[k]]
    if (verbose) message("Step ", k, ": ", length(parents), " parent models.")
    
    children_all <- list()
    child_keys   <- character()
    
    for (parent in parents) {
      
      parent_key <- model_key(parent)
      parent_aic <- aic_by_model[[parent_key]]
      
      unused <- setdiff(predictors, parent)
      if (length(unused) == 0) next
      
      # Generate child candidates by adding each unused variable
      child_info <- lapply(unused, function(v) {
        vars <- c(parent, v)
        key  <- model_key(vars)
        
        # compute AIC once
        if (is.null(aic_by_model[[key]])) {
          aic_by_model[[key]] <<- model_aic(vars)
        }
        
        list(vars = vars, key = key, AIC = aic_by_model[[key]])
      })
      
      # Evaluate children
      child_AICs <- sapply(child_info, `[[`, "AIC")
      best_child <- min(child_AICs)
      
      # Keep children: must (1) improve parent by eps and (2) be within delta
      keep_idx <- which(
        (parent_aic - child_AICs >= eps) &
          (child_AICs - best_child <= delta)
      )
      
      if (length(keep_idx) > 0) {
        kept <- child_info[keep_idx]
        children_all <- c(children_all, lapply(kept, `[[`, "vars"))
        child_keys   <- c(child_keys, sapply(kept, `[[`, "key"))
      }
    }
    
    # No children → stop early
    if (length(child_keys) == 0) {
      if (verbose) message("Stopping: no AIC-improving children at step ", k, ".")
      break
    }
    
    # Deduplicate
    uniq_index <- !duplicated(child_keys)
    frontier <- children_all[uniq_index]
    
    # Prune if too many
    if (length(frontier) > L) {
      AICs <- sapply(frontier, function(vars) aic_by_model[[model_key(vars)]])
      frontier <- frontier[order(AICs)[1:L]]
      if (verbose) message("Pruning to best ", L, " models.")
    }
    
    model_tree[[k + 1]] <- frontier
  }
  
  return(list(
    models_by_step = model_tree,
    aic_by_model   = aic_by_model,
    meta = list(
      family = family,
      K = K,
      eps = eps,
      delta = delta,
      L = L
    )
  ))
}

