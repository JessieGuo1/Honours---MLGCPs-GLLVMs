mlgcp_alt_glm = function(u, data, Z, X, N, max_it, occ, wts, 
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
  
  time_i = system.time({for (iter in 1:N){
    # Covariates and current spatial field
    if (one.quad == TRUE) {
      # Covariates and current spatial field
      X_pres_Zu = cbind(X_pres, Z_pres_u)
      X_quad_Zu = cbind(X_quad, Z_quad_u)
    }
    
    else {
      X_Zu = cbind(X, Zu)
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
      
      start_iter = if (all(is.finite(beta_lambda[j, ]))) {
        beta_lambda[j, ]
      } else {
        NULL
      }
      
      # Weighted Poisson ridge regression - DWPR response (pt/wt), weight = wt
      glm1 = glm.fit(y = occ_i/w_i, x = X_i, family = poisson(link = 'log'), weights = w_i, intercept = TRUE, start = start_iter, control = glm.control(maxit = max_it))
      
      # Store this species' intercept + beta + lambda (1 loading)
      beta_lambda[j, ] = t(c(glm1$coefficients)) 
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

    if (basis.function.update == TRUE) {
      # Basis function update for u
      K = nrow(u)
      D = ncol(u)
      offset_i = offset + Zl_joint %*% as.vector(u)
      for (i in 1:nrow(u)){
        idx = i + (0:(D - 1)) * K
        offset_i = offset_i - Zl_joint[, idx, drop = FALSE] %*% u[i, ] #
        start_iter_u = if (all(is.finite(u[i, ]))) {
          u[i, ]} else {
            NULL
          }        
        glm1 = glm.fit(x = Zl_joint[, idx, drop = FALSE], y = occ/wts, intercept = FALSE, family=poisson(link='log'), weights = wts, offset = as.numeric(offset_i), start = start_iter_u, control = glm.control(maxit = max_it)) 
        
        # Update shared-field coefficients
        u[i, ] = glm1$coefficients
        #u[i, ] = u[i, ] + new_u[i, ]
        #offset_i = offset_i + Zl_joint[, i, drop = FALSE]%*% new_u[i, ]
        offset_i = offset_i + Zl_joint[, idx, drop = FALSE]%*% u[i, ]    
      }
    }
    
    else {
      start_iter_u = if (all(is.finite(u))) {
        as.vector(u)} else {
          NULL
        }
      # Joint update for u
      glm1 = glm.fit(x = Zl_joint, y = occ/wts, intercept = FALSE, family=poisson(link='log'), weights = wts, offset = offset, start = start_iter_u, control = glm.control(maxit = max_it)) 
      
      # Update shared-field coefficients
      u <- matrix(coef(glm1), nrow = nrow(u), ncol = ncol(u))
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
    pred_mglm = mu
    
    # Weighted Poisson-process log-likelihood (DWPR/Berman-Turner approximation)
    loglik_new_mglm = sum(occ * eta - wts * exp(eta))
    cat("iter", iter, "loglik:", loglik_new_mglm, "\n")
    
    # Stop once relative change in log-likelihood is small
    if (is.finite(loglik) && abs(loglik_new_mglm - loglik)/ max(1, abs(loglik_new_mglm)) < 1e-6) break
    loglik = loglik_new_mglm}
  }) 
  return (list(pred_mglm = pred_mglm, eta = eta, loglik_new_mglm = loglik_new_mglm, time_mglm = as.numeric(time_i[3]), u = u, lambda = lambda, beta = beta, ite_converge = iter))
}
