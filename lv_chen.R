# Constrained JMLE fit
lv_chen = function(sim_counts, X, u, K, N){

  c = 5 * sqrt(K) # constraint radius
  max_it = 3
  
  # Storage for projected updates
  beta_proj = matrix(NA, nrow = ncol(sim_counts), ncol = 1 + ncol(X)) 
  lambda_proj = matrix(NA, nrow = ncol(sim_counts), ncol = ncol(u)) 
  u_proj = matrix(NA, nrow = nrow(sim_counts), ncol = ncol(u)) #56 x 1
  
  # Initialise lik
  loglik = -Inf
  new_X = cbind(1, X) #56 x 3
  
  # initialize beta, lambda from Poisson GLMs
  time_i = system.time({
    for (iter in 1:N) {
      
      # Column update, update beta, lambda
      X_u = cbind(new_X, u)
      
      suppressWarnings(
        for (j in 1:ncol(sim_counts)) {
          # Warm starts
          if (iter == 1) {
            warm_start_j = NULL
          }
          else {
            warm_start_j = as.numeric(beta_lambda[j, ])
          }        
          
          glm1 = try(glm.fit(y = sim_counts[, j], x = X_u, intercept = TRUE, family = poisson(link = "log"), control = glm.control(maxit = max_it), start = warm_start_j), silent = TRUE)
          
          if (inherits(glm1, "try-error") || any(!is.finite(glm1$coefficients))) {
            glm1 = glm.fit(y = sim_counts[, j], x = X_u, intercept = TRUE, family = poisson(link = "log"), start = NULL, control = glm.control(maxit = max_it))
          }
          
          new_beta_lambda = as.numeric(t(glm1$coefficients))
          
          # Extract beta, lambda
          beta_j = new_beta_lambda[1:(1 + ncol(X)), drop = FALSE] #1 x 3
          lambda_j = new_beta_lambda[(ncol(X) + 2):(ncol(X) + 1 + ncol(u)), drop = FALSE] #1 x 1
          
          # Constrain beta, lambda
          comb_j = c(beta_j, lambda_j) #4 x 1
          comb_proj = proc(comb_j, c)
          
          beta_proj[j, ] = comb_proj[1:(ncol(X) + 1), drop = FALSE]
          lambda_proj[j, ] = comb_proj[(ncol(X) + 2):(ncol(X) + 1 + ncol(u)), drop = FALSE]
        }
      )
      
      # Update parameters with constrained versions
      beta_lambda = cbind(beta_proj, lambda_proj)
      offset = new_X%*%t(beta_proj)
      
      # Row update
      suppressWarnings(
        for (i in 1:nrow(sim_counts)) {
          
          # Fit u with offset 
          glm1 = try(glm.fit(y = sim_counts[i, ], x = lambda_proj, intercept = FALSE, family = poisson(link = "log"), offset = offset[i, ], control = glm.control(maxit = max_it), start = as.numeric(u[i, ])), silent = TRUE)
          
          if (inherits(glm1, "try-error") || 
              any(!is.finite(glm1$coefficients))) {
            glm1 = glm.fit(y = sim_counts[i, ], x = lambda_proj, intercept = FALSE, family = poisson(link = "log"), offset = offset[i, ], control = glm.control(maxit = max_it), start = NULL)
          }
          
          new_u = as.matrix(t(glm1$coefficients))
          # Constrain u
          u_proj[i, ] = proc(new_u, sqrt(c^2 - 1))
        }
      )
      u = u_proj
      
      # Fitted means, log-lik
      eta = offset + u %*% t(lambda_proj) #56 x 200
      #      eta = pmin(pmax(eta, -10), 10)
      mu = exp(eta)
      pred_counts_i = mu
      
      loglik_new = sum(dpois(sim_counts, lambda = mu, log = TRUE))
      # Convergence check        
      if (is.finite(loglik) && is.finite(loglik_new)) {
        rel_change = abs(loglik_new - loglik)/ max(1, abs(loglik_new)) 
        
        cat(
          "iter =", iter,
          "loglik =", loglik,
          "beta = ", range(beta_proj),
          "lambda = ", range(lambda_proj),
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
  return(list(time = as.numeric(time_i[3]), ite_converge = iter, loglik_new = loglik_new, beta = beta_proj, lambda = lambda_proj, u = u, pred_counts = pred_counts_i, eta = eta))
}