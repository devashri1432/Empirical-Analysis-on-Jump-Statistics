# Load necessary libraries
library(forecast)
library(cluster)
library(stats)
library(SLBDD)
source("JUMPMETHODFUNCTIONS.R")

# Function generating DGP-4
dgp_4 <- function(n, df, cluster_size, discard = 20) {
  total_no_series <- sum(cluster_size)
  total_length <- n + discard
  
  data_matrix <- matrix(NA, nrow = n, ncol = total_no_series)
  
  for (i in 1:cluster_size[1]) {
    
    ts <- arima.sim(n = total_length, model = list(ar = c(1.3, -0.4)), innov = rt(total_length,df))
    data_matrix[, i] <- ts[(discard + 1):total_length]
  }
  
  for (i in (cluster_size[1]+1):(total_no_series-cluster_size[3])) {
    
    ts <- arima.sim(n = total_length, model = list(ar = c(0.8)), innov = rt(total_length,df))
    data_matrix[, i] <- ts[(discard + 1):total_length]
  }
  
  for (i in ((total_no_series-cluster_size[3])+1):total_no_series) {
    
    ts <- arima.sim(n = total_length, model = list(ma = c(-0.1, -0.42)), innov = rt(total_length,df))
    data_matrix[, i] <- ts[(discard + 1):total_length]
  }
  return(data_matrix)
}

num_series <- 33  
sample_size <- 300  
max_clusters <- 10 
num_iter <- 500
cluster_size_sim4 <- c(10,15,8)
df <- 10
count_jump <- 0
true_no_of_clust <- 3
optimal_no_of_clust <- 0

for (iteration in 1:num_iter) {
  set.seed(iteration)
  
  data_matrix <- dgp_4(sample_size, df, cluster_size_sim4)

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

print("Empirical Probability for Simulation-4 Jump Statistic:")
emp_prob <- count_jump / num_iter
print(emp_prob)
