check_stationarity <- function(ar_coeff) {
  roots <- polyroot(c(1, -ar_coeff))  # Find roots of the characteristic equation
  return(all(Mod(roots) > 1))  # Check if all roots lie outside the unit circle
}

generate_stationary_ar <- function(order) {
  while (TRUE) {
    ar_coeff <- runif(order, -1, 1)  # Generate random AR coefficients
    if (check_stationarity(ar_coeff)) {
      return(ar_coeff)
    }
  }
}

J_original_distribution <- function(X) {
  
  # if (is.null(X) || !is.matrix(X) || ncol(X) < 2) {
  #   stop("Input matrix X must be valid and have at least 2 columns.")
  # }
  
  k <- ncol(X)
  
  if (!is.numeric(k) || k <=0 ){
    stop("Error in generating J_0!")
  }
  pacf_features <- matrix(0, nrow = k, ncol = 5)
  for (i in 1:k) {
    pacf_result <- pacf(X[, i], plot = FALSE, lag.max = 5)
    pacf_features[i, ] <- pacf_result$acf[1:5]
  }
  
  dist_matrix <- as.matrix(dist(pacf_features), method = "euclidean")
  
  hc <- hclust(as.dist(dist_matrix), method = "complete")
  dendrogram_heights <- hc$height
  # print(dendrogram_heights)
  jumps <- diff(dendrogram_heights)
  # print(length(jumps))
  
  if (length(jumps) < (k - 2)) {
    jumps <- c(jumps, rep(0, (k - 2) - length(jumps)))  # Pad with zeros
  } else if (length(jumps) > (k - 2)) {
    jumps <- jumps[1:(k - 2)]  # Truncate to (k - 2)
  }
  
  J_0 <- jumps
  return(J_0)
}


J_reference_bootstrap_distribution <- function(X, B=100){
  
  T <- nrow(X)
  k <- ncol(X)
  
  l1_distances <- colSums(sapply(1:ncol(X), function(i) colSums(abs(X - X[, i]))))
  median_index <- which.min(l1_distances)
  x_ut <- X[, median_index] # Median Empirical Dynamic Quantile (Median EDQ)
  
  ar_model <- auto.arima(x_ut, max.q = 0, stationary = TRUE, seasonal = FALSE) # This generates only stationary data
  ar_coeff <- ar_model$coef # Extract AR coefficients
  
  # print(ar_coeff)
  
  residuals <- residuals(ar_model)
  ahat_ut <- residuals - mean(residuals, na.rm = TRUE)
  # print(ahat_ut)
  
  # print(centered_residuals)
  
  J_k <- matrix(0, nrow = k-2, ncol = B)
  X_b <- matrix(0, nrow = T, ncol = k)
  
  for (b in 1:B) {
    
    if (!check_stationarity(ar_coeff)) {
      cat("Non-stationary AR coefficients encountered at iteration", b, ". Replacing with stationary coefficients.\n")
      ar_coeff <- generate_stationary_ar(length(ar_coeff))
    }
    
    if (!is.numeric(k) || k <=0 ){
      stop("Error in generating J_k!")
    }
    
    pacf_features <- matrix(0, nrow = k, ncol = 5)
    for (j in 1:k) {
      # Sample T values with replacement from centered_residuals
      bootstrap_resid <- sample(ahat_ut, T, replace = TRUE)
      # Generate time series of length T using the AR coefficients
      X_b[, j] <- arima.sim(n = T, model = list(ar = ar_coeff), innov = bootstrap_resid)
    }
    
    
    # Step 2: Calculate ACF features from bootstrap samples
    for (j in 1:k) {
      # Calculate ACF for series starting from the 21st observation
      pacf_result <- pacf(X_b[, j], plot = FALSE, lag.max = 5)
      pacf_features[j, ] <- pacf_result$acf[1:5]
    }
    
    # Step 3: Calculate Euclidean distance and perform hierarchical clustering
    dist_matrix <- as.matrix(dist(pacf_features, method = "euclidean"))
    
    hc <- hclust(as.dist(dist_matrix), method = "complete")
    
    # Step 4: Calculate jumps from dendrogram heights and store in jump matrix
    dendrogram_heights <- hc$height  
    jumps <- diff(dendrogram_heights)
    
    # if (length(jumps) < (k - 2)) {
    #   jumps <- c(jumps, rep(0, (k - 2) - length(jumps)))  # Pad with zeros
    # } else if (length(jumps) > (k - 2)) {
    #   jumps <- jumps[1:(k - 2)]  # Truncate to (k - 2)
    # }
    
    # Store the adjusted jumps in J_k
    J_k[, b] <- jumps
  }
  return(J_k)
}

Jump_test <- function(J_0, J_k, alpha_values=c(0.01,0.025,0.05), B = 100){
  
  T_values <- numeric(length(alpha_values))
  # Calculate the reference distribution of critical values C(alpha1)
  critical_values <- numeric(length(alpha_values))
  
  # Calculate T(alpha) for each alpha in alpha_values
  for (j in 1:length(alpha_values)) {
    alpha <- alpha_values[j]
    # Calculate the 1-alpha quantile of the observed jumps in J for T(alpha)
    T_values[j] <- quantile(J_0, probs = 1 - alpha)
  }
  
  for (j in 1:length(alpha_values)) {
    alpha <- alpha_values[j]
    # Compute 1-alpha quantile for each bootstrap sample
    bootstrap_quantiles <- apply(J_k, 2, function(x) quantile(x, probs = 1 - alpha))
    # C(alpha1) is the 1-alpha1 quantile of these bootstrap quantiles
    critical_values[j] <- quantile(bootstrap_quantiles, probs = 1 - alpha, na.rm = TRUE)
  }
  
  # Decision rule: reject H0 if any T(alpha) > C(alpha1)
  reject_H0 <- any(T_values > critical_values, na.rm = TRUE)
  
  conclusion <- if (isTRUE(reject_H0)) {
    "Reject H0: There is evidence of more than one cluster (Ha: G > 1)"
  } else {
    "Fail to reject H0: Evidence supports a single cluster (H0: G = 1)"
  }
  
  # Return results
  test_results = list(
    T_values = T_values,
    critical_values = critical_values,
    conclusion = conclusion,
    reject_H0 = reject_H0
  )
  return(test_results)
}

calculate_test_statistic <- function(cluster_data, alpha_values = c(0.01, 0.025, 0.05), B = 100) {
  
  # if (is.null(cluster_data) || !is.matrix(cluster_data) || ncol(cluster_data) < 2) {
  #   print("Skipping cluster due to insufficient data.")
  #   return(FALSE)  # Treat as homogeneous if insufficient data
  # }
  
  J_0 <- J_original_distribution(cluster_data)  
  J_k <- J_reference_bootstrap_distribution(cluster_data, B = B)
  
  # Use perform_jump_test function to check homogeneity of the cluster
  test_results <- Jump_test(J_0, J_k, alpha_values = alpha_values, B = B)
  
  # Return whether to reject the null hypothesis of a single cluster
  return(test_results$reject_H0)
}

determine_optimal_no_of_clusters <- function(X, Gmax, alpha_values = c(0.01, 0.025, 0.05), B = 100){
  
  g_hat <- 2
  k <- ncol(X)
  
  if (!is.numeric(k) || k <=0 ){
    stop("Error in determining optimal number of clusters!")
  }
  
  pacf_features <- matrix(0, nrow = k, ncol = 5)
  for (i in 1:k) {
    pacf_result <- pacf(X[, i], plot = FALSE, lag.max = 5)
    pacf_features[i, ] <- pacf_result$acf[1:5]
  }
  
  dist_matrix <- as.matrix(dist(pacf_features, method = "euclidean"))
  
  repeat {
    # Step 3: Perform hierarchical clustering with g_hat clusters
    hc <- hclust(as.dist(dist_matrix), method = "complete")
    cluster_labels <- cutree(hc, k = g_hat)  # Cut tree to form g_hat clusters
    
    
    # Step 4: Check each cluster for homogeneity
    homogeneous <- TRUE
    for (j in 1:g_hat) {
      # Extract data for cluster j
      cluster_data <- X[, cluster_labels == j, drop = FALSE]
      
      # Skip clusters with insufficient data
      if (ncol(cluster_data) < 2) {
        next
      }
      
      if (calculate_test_statistic(cluster_data, alpha_values = alpha_values, B = B)) {
        homogeneous <- FALSE
        break
      }
    }
    
    # If all clusters are homogeneous, we've found the optimal number of clusters
    if (homogeneous || g_hat >= Gmax) {
      break
    }
    
    # Otherwise, increment g_hat and repeat the process
    g_hat <- g_hat + 1
  }
  
  return(g_hat)
}