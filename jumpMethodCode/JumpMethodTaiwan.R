check_stationarity <- function(ar_coeff) {
  roots <- polyroot(c(1, -ar_coeff))
  return(all(Mod(roots) > 1))
}

generate_stationary_ar <- function(order) {
  while (TRUE) {
    ar_coeff <- runif(order, -1, 1)
    if (check_stationarity(ar_coeff)) {
      return(ar_coeff)
    }
  }
}

# Function to compute observed set of jumps for original data
ObservedJumps <- function(X) {
  
  k <- ncol(X)
  n <- 6# No. of desired first sample ACF values. Modify as per the requirements of analysis 
  
  if (!is.numeric(k) || k <= 0){
    stop("ERROR: Observed jumps cannot be generated!")
  }
  
  acf_features_matrix <- matrix(NA, nrow = k, ncol = n)
  for (i in 1:k) {
    acf_result <- acf(X[, i], plot = FALSE, lag.max = n)
    acf_features_matrix[i, ] <- acf_result$acf[2:(n+1)] # Ignoring lag 0
  }
  
  dist_matrix <- as.matrix(dist(acf_features_matrix, method = "euclidean"))
  
  hc_result <- hclust(as.dist(dist_matrix), method = "complete")
  dend_heights  <- hc_result$height
  # print(dend_heights )

  J_0 <- diff(dend_heights)
  # print(length(jumps))
  
  return(J_0)
}

# Function to compute reference set of jumps from the null reference distribution obtained for original data
ReferenceJumps <- function(X, B=100){
  
  T <- nrow(X)
  k <- ncol(X)
  n <- 6
  
  l1_dist <- colSums(sapply(1:ncol(X), function(i) colSums(abs(X - X[, i]))))
  median_EDQ_index <- which.min(l1_dist)
  x_ut <- X[, median_EDQ_index] # Median EDQ of the original dataset
  
  adf_test <- suppressWarnings(tseries::adf.test(x_ut))
  
  if (adf_test$p.value < 0.05) {
    ar_model <- auto.arima(x_ut, max.q = 0, max.d = 0, stationary = TRUE, seasonal = FALSE)
  } else {
    ar_model <- auto.arima(x_ut, max.q = 0, max.d = 0, seasonal = FALSE)
  }
  
  ar_model <- auto.arima(x_ut, max.q = 0, stationary = TRUE, seasonal = FALSE)
  ar_coeff <- ar_model$coef
  # print(ar_coeff)
  
  residuals <- residuals(ar_model)
  ahat_ut <- residuals - mean(residuals, na.rm = TRUE)
  # print(ahat_ut)
  # print(centered_residuals)
  
  J_k <- matrix(NA, nrow = k-2, ncol = B)
  
  if (!is.numeric(k) || k <=0){
    stop("ERROR: Reference jumps cannot be generated!")
  }
  
  for (b in 1:B) {
    if (!check_stationarity(ar_coeff)) {
      cat("Non-stationary AR coefficients found at iteration", b, ". Replacing with stationary AR coefficients.\n")
      ar_coeff <- generate_stationary_ar(length(ar_coeff))
    }
    
    X_b <- matrix(NA, nrow = T, ncol = k)
    for (j in 1:k) {
      bootstrap_resid <- sample(ahat_ut, T, replace = TRUE)
      X_b[, j] <- arima.sim(n = T, model = list(ar = ar_coeff))
    }
    
    acf_features_matrix <- matrix(NA, nrow = k, ncol = n)
    for (j in 1:k) {
      acf_result <- acf(X_b[, j], plot = FALSE, lag.max = n)
      acf_features_matrix[j, ] <- acf_result$acf[2:((n+1))]
    }
    
    dist_matrix <- as.matrix(dist(acf_features_matrix, method = "euclidean"))
    
    hc_result <- hclust(as.dist(dist_matrix), method = "complete")
    
    dend_heights  <- hc_result$height  
    ref_jumps <- diff(dend_heights)
    
    J_k[, b] <- ref_jumps
  }
  return(J_k)
}

# Function to perform Jump test
JumpTest <- function(J_0, J_k, alpha_values=c(0.01,0.025,0.05), B = 100){

  T_values <- numeric(length(alpha_values))
  critical_values <- numeric(length(alpha_values))
  
  for (j in 1:length(alpha_values)) {
    
    alpha <- alpha_values[j]
    T_values[j] <- quantile(J_0, probs = 1 - alpha, na.rm = TRUE)
    boot_quantiles <- apply(J_k, 2, function(x) quantile(x, probs = 1 - alpha, na.rm = TRUE))
    critical_values[j] <- quantile(boot_quantiles, probs = 1 - alpha, na.rm = TRUE)
}
  
  reject_H0 <- any(T_values > critical_values, na.rm = TRUE)
  conclusion <- if (isTRUE(reject_H0)) {
    "Reject H0: There are multiple clusters (H1: G > 1)"
    } 
    else {
    "Fail to Reject H0: There is only single cluster (H0: G = 1)"
    }
  
  test_result_list = list(
    T_values = T_values,
    critical_values = critical_values,
    conclusion = conclusion,
    reject_H0 = reject_H0
  )
  return(test_result_list)
}

TestStatisticCalculationforClusterData <- function(clust_data, alpha_values = c(0.01, 0.025, 0.05), B = 100) {
  
  if (is.null(clust_data) || !is.matrix(clust_data) || ncol(clust_data) < 2) {
    return(list(reject_H0 = FALSE, reason = "Not enough data for testing"))
  }
  
  J_0 <- ObservedJumps(clust_data)  
  J_k <- ReferenceJumps(clust_data, B = B)
  
  test_result_list <- Jump_test(J_0, J_k, alpha_values = alpha_values, B = B)
  
  return(test_result_list$reject_H0)
}

OptimalNoOfClus <- function(X, G_max, alpha_values = c(0.01, 0.025, 0.05), B = 100){
  
  if (!is.numeric(G_max) || G_max <= 1) stop("G_max must be greater than 1!")
  if (!is.numeric(alpha_values) || any(alpha_values <= 0 | alpha_values >= 1)) stop("alpha_values must be between 0 and 1!")
  if (!is.matrix(X)) stop("X must be a matrix!")
  
  g_hat <- 2
  k <- ncol(X)
  n <- 6
  
  if (!is.numeric(k) || k <= 0){
    stop("Error in determining optimal number of clusters!")
  }
  
  acf_features_matrix <- matrix(NA, nrow = k, ncol = n)
  for (i in 1:k) {
    acf_result <- acf(X[, i], plot = FALSE, lag.max = n)
    acf_features_matrix[i, ] <- acf_result$acf[2:(n+1)]
  }
  
  dist_matrix <- as.matrix(dist(acf_features_matrix, method = "euclidean"))
  
  repeat {
    hc_result <- hclust(as.dist(dist_matrix), method = "complete")
    # Cut dendrogram tree to for g_hat clusters
    clust_labels <- cutree(hc_result, k = g_hat) 
    
    homo <- TRUE
    for (j in 1:g_hat) {
      # Extracting data for cluster j and performing jump test for every cluster j
      clust_data <- X[, clust_labels == j, drop = FALSE]
      test_result <- TestStatisticCalculationforClusterData(clust_data, alpha_values = alpha_values, B = B)
      
      if (isTRUE(test_result)) {
        homo <- FALSE
        break
      }
    }
    
    if (homo || g_hat >= G_max) {
      break
    }
    
    g_hat <- g_hat + 1
  }
  
  return(g_hat)
}