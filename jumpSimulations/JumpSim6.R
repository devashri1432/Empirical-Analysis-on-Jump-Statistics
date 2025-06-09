# Load necessary libraries
library(forecast)
library(cluster)
library(stats)
library(SLBDD)
source("JUMPMETHODFUNCTIONS.R")

# Function generating DGP-6
dgp_6 <- function(n, cluster_size, discard = 20) {
  data_matrix <- matrix(NA, nrow = n, ncol = cluster_size)
  total_length <- n + discard
  
  for (i in 1:cluster_size) {
    ts <- arima.sim(n = total_length, model = list(ar = c(1.3, -0.4)), innov = rnorm(total_length))
    data_matrix[, i] <- ts[(discard + 1):total_length]
  }
  return(data_matrix)
}

cluster_size_sim6 <- 200  
sample_size <- 300  
max_clusters <- 10
num_iter <- 500
true_no_of_clust <- 1
count_jump <- 0
optimal_no_of_clust <- 0

for (iteration in 1:num_iter) {
  set.seed(iteration)
  
  data_matrix <- dgp_6(sample_size, cluster_size_sim6)

  # acf_features <- matrix(NA, nrow = num_series, ncol = 5)
  # for (i in 1:num_series) {
  #   acf_result <- acf(data_matrix_discarded[, i], plot = FALSE, lag.max = 5)
  #   acf_features[i, ] <- acf_result$acf[2:6]
  # }
  # 
  # dist_matrix <- as.matrix(dist(acf_features, method = "euclidean"))
  
  J_0 <- J_original_distribution(data_matrix)
  
  J_k <- J_reference_bootstrap_distribution(data_matrix, B = 100)
  
  results <- Jump_test(J_0, J_k)
  
  if (results$conclusion == "Reject H0: There is evidence of more than one cluster (Ha: G > 1)"){
    optimal_no_of_clust <- determine_optimal_no_of_clusters(data_matrix, Gmax=max_clusters)
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

print("Empirical Probability for Simulation-1 Jump Statistic:")
emp_prob <- count_jump / num_iter
print(emp_prob)