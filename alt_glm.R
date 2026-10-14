# Alternating GLM fitting function
alt_glm = function(sim_counts, X, u, N){
  
  # Initial log-lik used for convergence check
  loglik = -Inf
  new_X = cbind(1, X) #56 x 3
  beta_lambda = matrix(NA, nrow = ncol(sim_counts), ncol = 1 + ncol(X) + ncol(u)) # Preset matrix of coefficients including betas and lambdas, nrow = no. of columns in the responses, ncol = beta0, beta1, beta2, lambda
  
  # Measure total elapsed fitting time
  time_i = system.time({
    
    # Alternate between species update and site updates until convergence or until N
    for (iter in 1:N){
      X_u = cbind(new_X, u)
      
      # Column update - estimate beta, lambda for each species, conditional on u
      #      suppressWarnings(
      for (j in 1:ncol(sim_counts)){
        if (iter == 1) {
          warm_start_j = NULL
        }
        else {
          warm_start_j = as.numeric(beta_lambda[j, ])
        }
        
        glm1 = try(
          glm.fit(y = sim_counts[, j], x = X_u, intercept = TRUE, family=poisson(link='log'), control = glm.control(maxit = max_it), start = warm_start_j), silent = TRUE) 
        # Fallback for when warm-start fit fails, retry using default starts
        if (inherits(glm1, "try-error") ||
            any(!is.finite(glm1$coefficients))) {
          glm1 = glm.fit(y = sim_counts[, j], x = X_u, intercept = TRUE, family=poisson(link='log'), start = NULL, control = glm.control(maxit = max_it))
        }
        
        beta_lambda[j, ] = t(glm1$coefficients) #1 x 4
      }
      #        )
      # Extract fixed-effect coefficients
      beta = beta_lambda[, 1:(ncol(X)+1), drop = FALSE] #200 x 3
      # Extract latent-variable loadings
      lambda = beta_lambda[, (ncol(X) + 2):(ncol(X) + 1 + ncol(u)), drop = FALSE] #200 x 1
      #      lambda = pmax(lambda, -2000)
      offset = new_X%*%t(beta)
      
      # Row update - estimate latent score u for each site, conditional on current beta, lambda 
      # Preset matrix for updated u
      new_u = matrix(NA, nrow = nrow(sim_counts), ncol = ncol(u)) #56 x 1
      
      #      suppressWarnings(
      
      for (i in 1:nrow(sim_counts)){
        
        glm1 = try(
          glm.fit(y = as.numeric(sim_counts[i, ]), x = lambda, intercept = FALSE, family=poisson(link='log'), offset = offset[i, ], control = glm.control(maxit = max_it), start = as.numeric(u[i, ])), silent = TRUE)  #invalid type without as.numeric
        # Fallback for when warm-started fit fails
        if (inherits(glm1, "try-error") ||
            any(!is.finite(glm1$coefficients))) {
          glm1 = glm.fit(y = sim_counts[i, ], x = lambda, intercept = FALSE, family=poisson(link='log'), offset = offset[i, ], control = glm.control(maxit = max_it), start = NULL)
        }
        
        new_u[i, ] = t(glm1$coefficients) 
      }
      #      )
      # Update u
      u = new_u #56 x 1
      # Fitted means, log-lik
      eta = offset + u%*%t(lambda) #56 x 200 + 56 x 200
      #      eta = pmax(eta, -2000) #clamping to limit the exp(eta)
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