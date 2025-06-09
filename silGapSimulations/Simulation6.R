# Load necessary libraries
library(forecast)
library(cluster)
library(stats)
library(SLBDD)

# Function generating DGP-1
dgp_6 <- function(n, cluster_size, discard = 20) {
  data_matrix <- matrix(NA, nrow = n, ncol = cluster_size)
  total_length <- n + discard
  
  for (i in 1:cluster_size) {
    ts <- arima.sim(n = total_length, model = list(ar = c(1.3, -0.4)), innov = rnorm(total_length))
    data_matrix[, i] <- ts[(discard + 1):total_length]
  }
  return(data_matrix)
}

num_series <- 200  
sample_size <- 300  
max_clusters <- 10 
num_iter <- 500
count_silhouette <- 0
count_gap <- 0
optimal_clusters_gap <- 0
optimal_clusters_sil <- 0
true_no_of_clust <- 1

silhouette_scores <- matrix(NA, nrow = num_iter, ncol = max_clusters)
gap_values <- matrix(NA, nrow = num_iter, ncol= max_clusters-1)

for (iteration in 1:num_iter) {
  set.seed(iteration)
  
  data_matrix <- dgp_6(sample_size, num_series)
  data_matrix_discarded <- data_matrix[-c(1:20), ]
  
  acf_features <- matrix(NA, nrow = num_series, ncol = 5)
  for (i in 1:num_series) {
    acf_result <- acf(data_matrix_discarded[, i], plot = FALSE, lag.max = 5)
    acf_features[i, ] <- acf_result$acf[2:6]
  }
  
  dist_matrix <- as.matrix(dist(acf_features, method = "euclidean"))
  
  # print(dim(dist_matrix))[1]
  
  # Silhouette Statistic Calculation
  silh_result <- silh.clus(nClus=max_clusters, distanceMatrix = as.dist(dist_matrix), method="complete")
  sil_result_coef <- silh_result$coef
  sil_result_row <- t(sil_result_coef$silIndex)
  # print(sil_result_row)
  
  optimal_clusters_sil <- silh_result$nClus
  
  # Gap Statistic Calculation
  
  sc1 <- hclust(as.dist(dist_matrix), method = "complete")
  memb <- cutree(sc1, k=1:max_clusters)
  gap_result <- gap.clus(DistanceMatrix = dist_matrix, Clusters = memb, B=100)
  # gap_result_ngroups <- gap_result$gap.values
  # gap_result_row <- t(gap_result_ngroups$gap)
  # print(gap_result_row)
  
  # gap_values[iteration,] <- gap_result_row
  
  optimal_clusters_gap <- gap_result$optim.k
  
  if (optimal_clusters_sil == true_no_of_clust){
    count_silhouette <- count_silhouette + 1
  }
  
  if (optimal_clusters_gap == true_no_of_clust){
    count_gap <- count_gap + 1
  }
  
  cat("Iteration:", iteration, "\n")
  cat("Optimal Silhouette Clusters:", optimal_clusters_sil, "\n")
  cat("Optimal Gap Clusters:", optimal_clusters_gap, "\n")
}

# print(data_matrix)
# print(acf_features)
# print(silhouette_scores)
# print(gap_values)

# Optimal Clusters using Silhouette Statistics

# avg_silhouette_scores <- colMeans(silhouette_scores, na.rm = TRUE)
# plot(1: max_clusters, avg_silhouette_scores, type = "b",
#        xlab = "No of Clusters", ylab = "Avg Silhouette Score",
#        main = "Simulation-1 Silhouette Statistic")

print("Empirical Probability for Silhouette Statistic:")
emp_prob_sil <- count_silhouette / num_iter
print(emp_prob_sil)

# Optimal Clusters using Gap Statistics

# avg_gap_values <- colMeans(gap_values, na.rm=TRUE)
# plot(1: (max_clusters-1), avg_gap_values, type = "b",
#       xlab = "No of Clusters", ylab = "Avg Gap",
#       main = "Simulation-1 Gap Statistic")

emp_prob_gap <- count_gap / num_iter
print("Empirical Probability for Gap Statistic:")
print(emp_prob_gap)


# acf(data_matrix[, 2], lag.max = 5, main = "ACF of First Time Series (First 5 Lags)")
# acf(data_matrix[, 3], lag.max = 5, main = "ACF of First Time Series (First 5 Lags)")
