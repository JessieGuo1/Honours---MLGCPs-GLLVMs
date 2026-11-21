mlgcp_alt_glmnet = function(u, data, Z, X, N, max_it, occ, wts, ridgePen, 
                      basis.function.update = FALSE, 
                      one.quad = FALSE, pres = NULL, quad = NULL, X_pres = NULL, X_quad = NULL, Z_pres = NULL, Z_quad = NULL, red.data = NULL){
  
  # Shared spatial field
  if (one.quad == TRUE) {
    pres_spid = as.numeric(str_extract(pres$spid, "\\d+"))
    # Shared spatial field
    Z_pres_u = Z_pres %*% u 
    Z_quad_u = Z_quad %*% u 
  }
  else {
    Zu = Z%*%u 
  }
  
  chosen_spid = unique(spid)
  n_species = length(chosen_spid)
  
  # Store species-specific intercepts, slopes and loadings
  beta_lambda = matrix(NA,  nrow = n_species, ncol = ncol(u) + ncol(X)) 
  
  # Row indices for each species
  if (one.quad == TRUE) {
    rows_by_species = split(seq_len(nrow(pres)), pres_spid)
    species_idx_pres = match(pres_spid, chosen_spid)
  }
  else {
    rows_by_species = split(seq_len(nrow(dat_rename)), spid)
  }
  
  # Species index for each row
  species_idx_row = match(spid, chosen_spid)
  
  # Outer (column-update + row-update) iterations
  loglik = -Inf
  
  # Species-specific warm starts
  warm_start = vector("list", n_species)
  
  # Penalise only the spatial loading
  penalty_factor = c(rep(0, ncol(X) - 1), rep(1, ncol(u)))
  
  time_i = system.time({for (iter in 1:N){
    # Covariates and current spatial field
    if (one.quad == TRUE) {
      # Covariates and current spatial field
      X_pres_Zu = cbind(X_pres[, -1], Z_pres_u)
      X_quad_Zu = cbind(X_quad[, -1], Z_quad_u)
    }
    
    else {
      X_Zu = cbind(X[, -1], Zu)
    }

    # Species update: for each species, fit fixed effects + basis-function loadings (lambda)
    
    for (j in seq_along(chosen_spid)){
      i = chosen_spid[j]
      idx = rows_by_species[[as.character(i)]]
      
      if (one.quad == TRUE){
        # Species-specific design matrix
        X_i = rbind(X_pres_Zu[idx, ], X_quad_Zu)
        # Species-specific responses
        occ_i = c(red.data$occ[idx], rep(0, nrow(quad)))
        # Species-specific weights
        w_i = c(red.data$wt[idx], quad$wt)
      }
      
      else {
      # Species-specific design matrix
      X_i = X_Zu[idx, ]
      occ_i = occ[idx]
      w_i = wts[idx]
      } 
      
      # Weighted Poisson ridge regression - DWPR response (pt/wt), weight = wt
      glmnet1 = glmnet:::glmnet.fit(y = occ_i/w_i, x = as.matrix(X_i), family = poisson(link = 'log'), weights = w_i, intercept = TRUE, warm = warm_start[[j]], alpha = 0, lambda = ridgePen, penalty.factor = penalty_factor, control = glmnet.control(mxitnr = max_it))
      
      # Store this species' intercept + beta + lambda (1 loading)
      beta_lambda[j, ] = t(c(glmnet1$a0, as.numeric(glmnet1$beta))) 
      
      # Update warm start
      warm_start[[j]] = list(a0 = glmnet1$a0, beta = as.numeric(glmnet1$beta))
    }
    
    # Split fitted coefficients into fixed effects (beta) and loadings (lambda)
    beta_lambda
    beta = beta_lambda[, 1:ncol(X)] 
    lambda = as.matrix(beta_lambda[, (ncol(X) + 1):ncol(beta_lambda)]) 
    
    # Expand loadings to observation rows
    lambda_rep = lambda[species_idx_row, , drop = FALSE]
    
    # Fixed-effect linear predictor
    offset = rowSums(X * beta[species_idx_row, , drop = FALSE])
    
    # Shared-field design matrix
    Zl_joint = do.call(cbind, lapply(seq_len(ncol(lambda_rep)),
             function(d) Z * lambda_rep[, d])
    )
    
    attr(Zl_joint, "xm") = rep(0, ncol(Zl_joint))
    attr(Zl_joint, "xs") = rep(1, ncol(Zl_joint))
    
    if (basis.function.update == TRUE) {
      # Basis function update for u
      K = nrow(u)
      D = ncol(u)
      offset_i = offset + Zl_joint %*% as.vector(u)
      for (i in 1:nrow(u)){
        idx = i + (0:(D - 1)) * K
        offset_i = offset_i - Zl_joint[, idx, drop = FALSE] %*% u[i, ] #
        
        glmnet1 = glmnet:::glmnet.fit(x = Zl_joint[, idx, drop = FALSE], y = occ/wts, intercept = FALSE, family=poisson(link='log'), weights = wts, offset = as.numeric(offset_i), alpha = 0, lambda = ridgePen, warm = list(a0 = 0, beta = u[i, ]), control = glmnet.control(mxitnr = max_it)) 
        
        # Update shared-field coefficients
        u[i, ] = as.numeric(glmnet1$beta)
        #u[i, ] = u[i, ] + new_u[i, ]
        #offset_i = offset_i + Zl_joint[, i, drop = FALSE]%*% new_u[i, ]
        offset_i = offset_i + Zl_joint[, idx, drop = FALSE]%*% u[i, ]    
      }
    }
    
    else {
      # Joint update for u
      glmnet1 = glmnet:::glmnet.fit(x = Zl_joint, y = occ/wts, intercept = FALSE, family=poisson(link='log'), weights = wts, offset = as.numeric(offset), alpha = 0, lambda = ridgePen, warm = list(a0 = 0, beta = as.vector(u)), control = glmnet.control(mxitnr = max_it)) 
      
      # Update shared-field coefficients
      u <- matrix(glmnet1$beta, nrow = nrow(u), ncol = ncol(u))
    }
    
    print(u)
    
    # Recompute shared spatial field with updated u, for next iteration's column update
    if (one.quad == TRUE){
      Z_pres_u = Z_pres %*% u 
      Z_quad_u = Z_quad %*% u
    }
    
    Zu = Z%*%u
    
    # Full linear predictor
    eta = offset + rowSums(Zu * lambda_rep)
    mu = exp(eta)
    pred_mglmnet = mu
    
    # Weighted Poisson-process log-likelihood (DWPR/Berman-Turner approximation)
    loglik_new_mglmnet = sum(occ * log(mu) - wts * mu)
    cat("iter", iter, "loglik:", loglik_new_mglmnet, "\n")
    
    # Stop once relative change in log-likelihood is small
    if (is.finite(loglik) && abs(loglik_new_mglmnet - loglik)/ max(1, abs(loglik_new_mglmnet)) < 1e-6) break
    loglik = loglik_new_mglmnet}
  }) 
  return (list(pred_mglmnet = pred_mglmnet, eta = eta, loglik_new_mglmnet = loglik_new_mglmnet, time_mglmnet = as.numeric(time_i[3]), u = u, lambda = lambda, beta = beta, ite_converge = iter))
}
