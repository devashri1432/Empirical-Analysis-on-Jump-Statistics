# Load necessary libraries
library(forecast)
library(cluster)
library(stats)
library(SLBDD)
source("JUMPMETHODFUNCTIONS.R")

# Function generating DGP-10
dgp_10 <- function(n, cluster_size) {
  total_no_series <- sum(cluster_size)
  data_matrix <- matrix(NA, nrow = n, ncol = total_no_series)
  
  for (i in 1:cluster_size[1]) {
    
    ts <- arima.sim(n = n, model = list(ar = c(1.4,-0.48)), innov = rnorm(n))
    data_matrix[, i] <- ts
  }
  
  for (i in (cluster_size[1]+1):(total_no_series-(cluster_size[3]+cluster_size[4]))) {
    
    ts <- arima.sim(n = n, model = list(ar = c(0.8)), innov = rnorm(n))
    data_matrix[, i] <- ts
  }
  
  for (i in ((total_no_series-(cluster_size[3]+cluster_size[4]))+1):(total_no_series-cluster_size[4])) {
    
    ts <- arima.sim(n = n, model = list(ma = c(0.1, -0.42)), innov = rnorm(n))
    data_matrix[, i] <- ts
  }
  
  for (i in ((total_no_series-cluster_size[4])+1):total_no_series) {
    
    ts <- arima.sim(n = n, model = list(ma = c(0.6)), innov = rnorm(n))
    data_matrix[, i] <- ts
  }
  return(data_matrix)
}

num_series <- 210  
sample_size <- 500  
max_clusters <- 10 
num_iter <- 500
cluster_size_sim10 <- c(40,70,60,40)
count_jump <- 0
true_no_of_clust <- 4
optimal_no_of_clust <- 0

for (iteration in 1:num_iter) {
  set.seed(iteration)
  
  data_matrix <- dgp_10(sample_size, cluster_size_sim10)
  data_matrix_discarded <- data_matrix[-c(1:20), ]
  
  # acf_features <- matrix(NA, nrow = num_series, ncol = 5)
  # for (i in 1:num_series) {
  #   acf_result <- acf(data_matrix_discarded[, i], plot = FALSE, lag.max = 5)
  #   acf_features[i, ] <- acf_result$acf[2:6]
  # }
  # 
  # dist_matrix <- as.matrix(dist(acf_features, method = "euclidean"))
  
  J_0 <- J_original_distribution(data_matrix_discarded)
  
  J_k <- J_reference_bootstrap_distribution(data_matrix_discarded, B = 100)
  
  results <- Jump_test(J_0, J_k)
  
  if (results$conclusion == "Reject H0: There is evidence of more than one cluster (Ha: G > 1)"){
    optimal_no_of_clust <- determine_optimal_no_of_clusters(data_matrix_discarded, Gmax=max_clusters)
  }
  else{
    optimal_no_of_clust <- 1
  }
  
  if (optimal_no_of_clust == true_no_of_clust){
    count_jump <- count_jump + 1
  }
  
  cat("Iteration:", iteration, "\n")
  cat("Test Statistic Values:", results$T_values, "\n")
  cat("Critical Values:", results$critical_values, "\n")
  cat("Conclusion:", results$conclusion, "\n")
  cat("Optimal Number of Clusters:", optimal_no_of_clust, "\n")
}

print("Empirical Probability for Simulation-10 Jump Statistic:")
emp_prob <- count_jump / num_iter
print(emp_prob)