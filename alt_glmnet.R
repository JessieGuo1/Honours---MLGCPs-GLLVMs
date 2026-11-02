# Alternating penalised Poisson GLM fit
alt_glmnet = function(sim_counts, X, u, N, pen_strength){
  
  loglik = -Inf
  new_X = cbind(1, X) #56 x 3
  max_it = 3
  
  # Preset matrix for beta, lambda for each species
  beta_lambda = matrix(NA, nrow = ncol(sim_counts), ncol = 1 + ncol(X) + ncol(u)) 
  # Warm starts
  warm_start = vector("list", ncol(sim_counts))
  warm_start_u = vector("list", nrow(sim_counts))
  
  time_i = system.time({
    for (iter in 1:N){
      
      X_u = cbind(X, u)
      # Column update, beta, lambda, beta unpenalised, lambda penalised
      #        suppressWarnings(
      for (j in 1:ncol(sim_counts)){
        
        if (j == 1 && ncol(u) > 1){
          X_u_j = X_u[, -ncol(X_u), drop = FALSE]
          penalty_j = c(rep(0, ncol(X)), rep(1, ncol(u) - 1))
        } else {
          X_u_j = X_u
          penalty_j = c(rep(0, ncol(X)), rep(1, ncol(u)))
        }
        
        glmnet1 = try(glmnet:::glmnet.fit(y = sim_counts[, j], x = X_u_j, weights = rep(1, nrow(sim_counts)), intercept = TRUE, family=poisson(link='log'), warm = warm_start[[j]], alpha = 0, lambda = pen_strength, penalty.factor = penalty_j, control = glmnet.control(mxitnr = max_it)), silent = TRUE)
        if (inherits(glmnet1, "try-error") || any(!is.finite(c(as.numeric(glmnet1$a0), as.numeric(glmnet1$beta))))) {
          glmnet1 = glmnet:::glmnet.fit(y = sim_counts[, j], x = X_u_j, weights = rep(1, nrow(sim_counts)), intercept = TRUE, family=poisson(link='log'), warm = NULL, alpha = 0, lambda = pen_strength, penalty.factor = penalty_j, control = glmnet.control(mxitnr = max_it))
        }
        
        # Store coefficients
        if (j == 1 && ncol(u) > 1) {
          beta_lambda[j, 1:(ncol(beta_lambda) - 1)] = t(c(glmnet1$a0, as.numeric(glmnet1$beta))) #1 x ncol(1+X_u)  as.numeric and c() because it's stored as sparse matrix
          beta_lambda[j, ncol(beta_lambda)] = 0
        } else {
          beta_lambda[j, ] = t(c(glmnet1$a0, as.numeric(glmnet1$beta)))
        }
        warm_start[[j]] = list(a0 = as.numeric(glmnet1$a0), beta = as.numeric(glmnet1$beta))
      }
      
      beta = beta_lambda[, 1:(ncol(X)+1), drop = FALSE] #200 x 3
      lambda = beta_lambda[, (ncol(X) + 2):(ncol(X) + 1 + ncol(u)), drop = FALSE] #200 x 1
      
      offset = new_X%*%t(beta)
      # Preset matrix for u
      new_u = matrix(NA, nrow = nrow(sim_counts), ncol = ncol(u)) #56 x 1
      
      # Row update, u
      for (i in 1:nrow(sim_counts)){
        glmnet1 = try(glmnet:::glmnet.fit(y = as.numeric(sim_counts[i, ]), x = lambda, weights = rep(1, ncol(sim_counts)), intercept = FALSE, family=poisson(link='log'), offset = offset[i, ], alpha = 0, lambda = pen_strength, warm = warm_start_u[[i]], control = glmnet.control(mxitnr = max_it)), silent = TRUE) #invalid type without as.numeric
        if (inherits(glmnet1, "try-error") ||
            any(!is.finite(c(as.numeric(glmnet1$a0), as.numeric(glmnet1$beta))))) {
          glmnet1 = glmnet:::glmnet.fit(y = as.numeric(sim_counts[i, ]), x = lambda, weights = rep(1, ncol(sim_counts)), intercept = FALSE, family=poisson(link='log'), offset = offset[i, ], alpha = 0, lambda = pen_strength, warm = NULL, control = glmnet.control(mxitnr = max_it))
        }
        
        # Store updated u
        new_u[i, ] = t(as.numeric(glmnet1$beta))
        warm_start_u[[i]] = list(a0 = 0, beta = as.numeric(glmnet1$beta))
      }
      # Update u
      u = new_u #56 x 1
      
      # Fitted means, log-lik
      eta = offset + u%*%t(lambda) #56 x 200 + 56 x 200
      mu = exp(eta)
      pred_counts_i = mu
      
      loglik_new = sum(dpois(as.matrix(sim_counts), lambda = mu, log = TRUE))
      # Convergence check        
      if (is.finite(loglik) && is.finite(loglik_new)) {
        rel_change = abs(loglik_new - loglik)/ max(1, abs(loglik_new))
        
        cat(
          "iter =", iter,
          "loglik =", loglik,
          "beta = ", range(beta),
          "lambda = ", range(lambda),
          "u = ", range(u),
          "eta = ", range(eta),
          "loglik_new =", loglik_new,
          "rel_change =", rel_change,
          "\n")
        if (rel_change < 1e-6) break
      }
      loglik = loglik_new
      if (is.finite(loglik) == FALSE) {stop("Loglikelihood is infinite")}
    }
  }) 
  
  return(list(time = as.numeric(time_i[3]), ite_converge = iter, loglik_new = loglik_new, beta = beta, lambda = lambda, u = u, pred_counts = pred_counts_i, eta = eta))
}
