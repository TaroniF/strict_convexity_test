# ---------------------------------------------------
#  Kernels and derivatives
# ---------------------------------------------------

# Epanechnikov kernel function
epanechnikov_kernel <- function(u) {
  out <- numeric(length(u))
  ind <- abs(u) <= 1
  out[ind] <- 0.75 * (1 - u[ind]^2)
  out
}

# Derivative of the Epanechnikov kernel
epanechnikov_derivative <- function(u) {
  out <- numeric(length(u))
  ind <- abs(u) <= 1
  out[ind] <- -1.5 * u[ind]
  out
}

# Integral of Epanechnikov kernel (assuming |u|<=1)
epanechnikov_integral <- function(u) {
  0.5 - 0.75 * (u - u^3/3)
}

# ---------------------------------------------------
#  Bandwidth selection via balanced k‑fold CV
# ---------------------------------------------------

# Min bandwidth s.t. density estimator is positive everywhere
min_positive_bandwidth <- function(x, tol = .Machine$double.eps) {
  # x   : numeric vector of observations
  # tol : a tiny positive offset to ensure strict positivity
  if (!is.numeric(x) || length(x) < 2) {
    stop("`x` must be a numeric vector of length >= 2")
  }
  # 1. sort
  xs  <- sort(x)
  # 2. compute adjacent gaps
  gaps <- diff(xs)
  # 3. half of the largest gap
  h0   <- max(gaps) / 2
  # 4. bump by a tiny amount so the kernel at exactly u=1 isn't zero
  max(h0, xs[1], 1-xs[length(x)]) + tol
}

# Bandwidth selection via k‑fold CV using L2-error
optimal_bandwidth_K_folds <- function(x, y, kernel,
                                      lower = min_positive_bandwidth(x),
                                      upper = 1,
                                      k_folds = 5) {
  n <- length(x)
  fold_ids <- rep(seq_len(k_folds), length.out = n)
  fold_ids <- sample(fold_ids)
  
  cv_error <- function(h) {
    if (h <= 0) return(Inf)
    err <- 0
    for (k in seq_len(k_folds)) {
      test_idx  <- which(fold_ids == k)
      train_idx <- setdiff(seq_len(n), test_idx)
      m_test    <- length(test_idx)
      m_train   <- length(train_idx)
      
      # Nadaraya Watson estimator fit
      # distance matrix and kernel weights (m_test × m_train)
      K_vals <- kernel(outer(x[test_idx], x[train_idx], "-") / h)
      W      <- matrix(K_vals, nrow = m_test, ncol = m_train,
                       dimnames = list(NULL, NULL))
      
      denom <- rowSums(W)
      
      if (any(denom == 0)) return(Inf)
      
      numer <- W %*% matrix(y[train_idx], ncol = 1)
      preds <- as.numeric(numer) / denom
      
      if (anyNA(preds)) return(Inf)
      
      
      err   <- err + sum((y[test_idx] - preds)^2)
    }
    err
  }
  optimize(cv_error, lower = lower, upper = upper)$minimum
}

# ---------------------------------------------------
#  Nadaraya Watson Matrix weights for derivative estimator
# ---------------------------------------------------

nadaraya_watson_derivative_weights <- function(W, f_hat, W_prime, h) {
  t1 <- sweep(W_prime, 1, f_hat, "/")
  t2 <- sweep(W, 1, rowSums(W_prime) / f_hat^2, "*")
  (t1 - t2) / h
}

# ---------------------------------------------------
#  Birke–Dette inverse‑derivative estimator
# ---------------------------------------------------

birke_dette_inverse_derivative_estimator <- function(h, x_grid, fun_grid) {
  M <- length(fun_grid)
  
  # 1) sort
  ord <- order(fun_grid, method = "radix")
  fgs <- fun_grid[ord]
  
  # 2) term1: proportion of points < fgs[i] - h
  term1_s <- findInterval(fgs - h, fgs, left.open = TRUE) / M
  
  # 3) precompute window bounds for term2
  lo <- findInterval(fgs - h, fgs, left.open = TRUE) + 1
  hi <- findInterval(fgs + h, fgs, left.open = FALSE)
  
  # 4) inline kernel‐integral sums
  inv_h  <- 1 / h
  inv_h3 <- inv_h^3
  term2_s <- numeric(M)
  
  for (i in seq_len(M)) {
    li <- lo[i]; hi_i <- hi[i]
    if (li <= hi_i) {
      # differences v = fgs[j] - fgs[i]
      v <- fgs[li:hi_i] - fgs[i]
      n <- hi_i - li + 1
      
      # sum of epanechnikov_integral(u) over j:
      #   sum[0.5 - 0.75*u + 0.25*u^3] 
      # = 0.5*n -0.75*sum(u) +0.25*sum(u^3)
      # with u = v/h
      sum_v  <- sum(v)
      sum_v3 <- sum(v * v * v)
      term2_s[i] <- (0.5 * n
                     - 0.75 * (sum_v * inv_h)
                     + 0.25 * (sum_v3 * inv_h3)
      ) / M
    }
  }
  
  # 5) combine & build approxfun on the sorted grid
  phi_s <- term1_s + term2_s
  approxfun(fgs, phi_s, rule = 2, ties = "ordered")
}


# ---------------------------------------------------
#  Generalised inverse helper
# ---------------------------------------------------

geninv_interp <- function(x, y, u) {
  # Ensure inputs are numeric and of equal length
  stopifnot(is.numeric(x), is.numeric(y), length(x) == length(y))
  
  # Remove missing pairs (keeps indices aligned)
  keep <- !(is.na(x) | is.na(y))
  x <- x[keep]; y <- y[keep]
  stopifnot(length(x) > 1)   # need at least two points for interpolation
  
  # Sort by y (required by findInterval)
  ord <- order(y, method = "radix")
  x <- x[ord]; y <- y[ord]
  
  # Replace ties in y by tiny jitter to enforce strict monotonicity
  dup <- duplicated(y)
  if (any(dup)) {
    eps <- .Machine$double.eps * max(abs(y)) + 1e-12
    y[dup] <- y[dup] + cumsum(dup) * eps
  }
  
  # Now y is strictly increasing
  k <- findInterval(u, y)
  res <- numeric(length(u))
  inside <- which(k > 0 & k < length(y))
  if (length(inside)) {
    j <- k[inside]
    t <- (u[inside] - y[j]) / (y[j + 1] - y[j])
    res[inside] <- x[j] + t * (x[j + 1] - x[j])
  }
  res[u <= y[1]] <- x[1]
  res[u >= y[length(y)]] <- x[length(x)]
  res
}

# ---------------------------------------------------
#  Birke–Dette best L2 estimator
# ---------------------------------------------------

# Returns Estiamtor with derivative bd_derivative which minimizes L2-distance to m_unconstrained
birke_dette_average_estimator <- function(
    bd_derivative,
    m_unconstrained,
    x_grid
) {
  M = length(x_grid)
  dz = 1/M
  
  #Approximate the cumulative integral
  I = dz*cumsum(bd_derivative)
  
  # Create a matrix of integral differences:
  # For each u (row) and x (column), compute I(x) - I(u)
  I_diff <- outer(I, I, FUN = function(u, x) x - u)
  
  # Create a matrix where each row i contains unconstrained_estimator(x_grid[i])
  unconstrained_estimator_values_mat <- matrix(m_unconstrained, M, M, byrow = FALSE)
  
  # Now, m_c(x, u) = m(u) + (I(x) - I(u))
  m_c_mat <- unconstrained_estimator_values_mat + I_diff
  
  # Use the theoretical result of theorem 2
  m_cstar_values <- colMeans(m_c_mat)
  
  # Make a function that interpolates or approximates \hat{m}_C^*(x) on [0,1].
  # We can do a simple linear interpolation via 'approx'.
  m_cstar_function <- function(x) {
    approx(x_grid, m_cstar_values, xout = x, rule=2)$y
  }
  
  list(m_c = m_cstar_function)
}

# ---------------------------------------------------
#  Test statistic
# ---------------------------------------------------

test_statistic <- function(phi, mprime, x_grid) {
  mean((phi(mprime) - x_grid)^2)
}

# ---------------------------------------------------
#  Main convexity test
# ---------------------------------------------------

convexity_test <- function(x, y, h_r = NULL, B_out = 100, B_in = 100, alpha = 0.05,
                           h_d = NULL, M = 100) {
  # (x,y) points
  # h_r nadaraya-watson bandwidth, cross validation is performed is NULL
  # B_out, B_in bootstrap repetition for double bootstrap
  # alpha test level
  # h_d bandwidth for Birke - Dette inverse function estimator, if NULL h_d = h_r
  # M discretization points on [0,1]
  
  
  n <- length(x)
  K   <- epanechnikov_kernel
  K_p <- epanechnikov_derivative
  
  # Discretization of [0,1] to compute estimators
  x_grid  <- (seq_len(M) - 0.5) / M
  diffmat <- outer(x_grid, x, "-")
  
  # Set h_r if not provided through cross-validation
  if (is.null(h_r)) h_r <- optimal_bandwidth_K_folds(x, y, K)
  
  # Kernel matrices
  W_k_vals  <- K(diffmat / h_r)
  W_k       <- matrix(W_k_vals, nrow = M, ncol = n)
  
  # Unnormalized density estimator
  f_hat     <- rowSums(W_k)
  
  # Check to avoid points where the estimated density is zero
  zero_row <- f_hat == 0 | is.na(f_hat)
  if (any(zero_row)) {
    keep    <- !zero_row
    W_k     <- W_k[keep, , drop = FALSE]
    f_hat   <- f_hat[keep]
    x_grid  <- x_grid[keep]
    diffmat <- diffmat[keep, , drop = FALSE]
  }
  
  # Matrices for Nadaraya Watson estimators of regression function and its derivative
  W         <- sweep(W_k, 1, f_hat, "/")
  W_kp_vals <- K_p(diffmat / h_r)
  W_k_p     <- matrix(W_kp_vals, nrow = M, ncol = n)
  W_p       <- nadaraya_watson_derivative_weights(W_k, f_hat, W_k_p, h_r)
  
  # Nadaraya Watson estimators
  m_hat      <- as.numeric(W %*% matrix(y, ncol = 1))
  m_primehat <- as.numeric(W_p %*% matrix(y, ncol = 1))
  
  # Set h_d if not provided to h_r
  if (is.null(h_d)) h_d <- h_r
  
  # Constrained estiamate of regression function and its derivative
  bd_inv_fun <- birke_dette_inverse_derivative_estimator(h_d, x_grid, m_primehat)
  phi_vec    <- bd_inv_fun(m_primehat)
  m_prime_iso<- geninv_interp(m_primehat, phi_vec, x_grid)
  bd_fit     <- birke_dette_average_estimator(m_prime_iso, m_hat, x_grid)
  m_c_hat    <- bd_fit$m_c(x)
  
  # Test statistics computation
  Tn <- test_statistic(bd_inv_fun, m_primehat, x_grid)
  
  
  
  ### Double bootstrap to extract 
  
  # Residuals
  res <- y - m_c_hat
  # Center residuals
  res <- res - mean(res)
  
  p1 <- (sqrt(5) + 1) / (2 * sqrt(5))
  
  Critical_vals <- numeric(B_out)
  Tn_boot_out <- numeric(B_out)
  p_value_in <- numeric(B_out)
  
  
  for (b in seq_len(B_out)) {
    # Wild bootstraps multipliers
    omega  <- ifelse(runif(n) <= p1, (1 - sqrt(5))/2, (1 + sqrt(5))/2)
    
    # Bootstrapped values
    y_star <- m_c_hat + omega * res
    
    # Refit estimators using (x, y_star)
    m_star <- as.numeric(W %*% matrix(y_star, ncol = 1))
    m_prime_star <- as.numeric(W_p %*% matrix(y_star, ncol = 1))
    bd_inv_fun_star  <- birke_dette_inverse_derivative_estimator(h_d, x_grid, m_prime_star)
    phi_star_vec    <- bd_inv_fun_star(m_prime_star)
    m_prime_star_iso<- geninv_interp(m_prime_star, phi_star_vec, x_grid)
    bd_star_fit     <- birke_dette_average_estimator(m_prime_star_iso, m_star, x_grid)
    m_c_star_hat    <- bd_star_fit$m_c(x)
    
    Tn_boot_out[b] <- test_statistic(bd_inv_fun_star, m_prime_star, x_grid)
    
    Tn_boot_in <- numeric(B_in)
    
    # Second order residuals
    res_star <- y_star - m_c_star_hat
    # Center residuals
    res_star <- res_star - mean(res_star)
    
    
    for (i in seq_len(B_in)){
      # Inner bootstrap
      omega_star  <- ifelse(runif(n) <= p1, (1 - sqrt(5))/2, (1 + sqrt(5))/2)
      
      y_star_star <- m_c_star_hat + omega_star * res_star
      
      m_prime_star_star <- as.numeric(W_p %*% matrix(y_star_star, ncol = 1))
      bd_inv_fun_star_star  <- birke_dette_inverse_derivative_estimator(h_d, x_grid, m_prime_star_star)
      
      Tn_boot_in[i]   <- test_statistic(bd_inv_fun_star_star, m_prime_star_star, x_grid)
    }
    
    # 1-alpha critical value of inner bootstrap
    Critical_vals[b] <- quantile(Tn_boot_in, 1 - alpha, names = FALSE)
    p_value_in[b] <- (1+sum(Tn_boot_out[b] <= Tn_boot_in))/(1+B_in)
  }
  
  # Decision and outputs
  p_value_observed <- (1+sum(Tn <= Tn_boot_out))/(B_out+1)
  
  
  # Decision based on critical region
  decision_critical_region_median <- as.integer(Tn >= median(Critical_vals))
  decision_critical_region_mean <- as.integer(Tn >= mean(Critical_vals))
  decision_critical_region_quantile <- as.integer(Tn >= quantile(Critical_vals, 1 - alpha, names = FALSE))
  
  
  # Decision based of p-value distribution
  p_val    <- (1+sum(p_value_in <= p_value_observed))/(1+B_out)
  decision_leq <- as.integer(p_val <= alpha)
  decision_le  <- as.integer(p_val < alpha)
  
  # Continuity correction

  p_val_corrected <- (1+sum(p_value_in < p_value_observed) + 0.5 * sum(p_value_in == p_value_observed))/(1+B_out)
  decision_leq_corrected <- as.integer(p_val_corrected <= alpha)
  decision_le_corrected  <- as.integer(p_val_corrected < alpha)
  
  
  
  structure(list(Tn = Tn, h_r = h_r, h_d = h_d,
                 Critical_vals = Critical_vals,
                 p_value_in = p_value_in,
                 Tn_boot_out = Tn_boot_out,
                 p_value_observed = p_value_observed,
                 p_value = p_val, 
                 p_val_corrected = p_val_corrected,
                 decision_leq = decision_leq,
                 decision_le = decision_le,
                 decision_leq_corrected = decision_leq_corrected,
                 decision_le_corrected = decision_le_corrected,
                 
                 decision_critical_region_median = decision_critical_region_median,
                 decision_critical_region_mean = decision_critical_region_mean,
                 decision_critical_region_quantile = decision_critical_region_quantile,
                 
                 x_grid = x_grid,
                 unconstrained_estimator = m_hat,
                 birke_dette_estimator   = bd_fit$m_c(x_grid),
                 unconstrained_derivative_estimator = m_primehat,
                 constrained_derivative_estimator = m_prime_iso,
            
                 inv_derivative_function = bd_inv_fun),
            
            class = "convexity test")
}





