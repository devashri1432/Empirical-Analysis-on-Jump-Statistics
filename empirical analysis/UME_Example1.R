# Required libraries
# install.packages("dendextend")
# install.packages("tidyverse")

library(forecast)
library(cluster)
library(stats)
library(SLBDD)
library(dendextend)
library(tidyverse)
source("JUMPMETHODFUNCTIONS2.R")

data("UMEdata20002018")

data_matrix <- as.matrix(UMEdata20002018)
col_names <- colnames(data_matrix)

data_matrix <- as.matrix(data_matrix)

growth_rates <- diff(log(data_matrix))

max_clusters <- 3

# Plotting the time series of quaterly growth rates of domestic product
years <- seq(2000, by = 0.25, length.out = nrow(growth_rates))
years <- years[-1]  # Remove the first element to match the length of growth_rates


custom_colors <- c("blue", "red", "green", "purple", "orange", "brown", "pink", "cyan", "black", "yellow")
col_vector <- rep(custom_colors, length.out = ncol(growth_rates))


matplot(years, growth_rates, type = "l", lty = 1, col = col_vector,
        xlab = "Year", ylab = "Rate", main = "")


years_original_data <- seq(2000, by = 0.25, length.out = nrow(data_matrix))
matplot(years_original_data, data_matrix, type = "l", lty = 1, col = col_vector,
        xlab = "Year", ylab = "Rate", main = "")

# Add a legend for clarity
# legend("topright", legend = paste("Series", 1:ncol(growth_rates)), 
#        col = rainbow(ncol(growth_rates)), lty = 1, cex = 0.6)

num_series <- ncol(growth_rates)

pacf_features <- matrix(NA, nrow = num_series, ncol = 5)
for (i in 1: num_series) {
  pacf_result <- pacf(growth_rates[, i], plot = FALSE, lag.max = 5)
  pacf_features[i, ] <- pacf_result$acf[1:5]
}

# print(pacf_features)

# Computing the similarity matrix by calculating the euclidean distance between each series
dist_matrix <- as.matrix(dist(pacf_features, method = "euclidean"))

# Silhouette Statistic Calculation

silh_result <- silh.clus(nClus = max_clusters, distanceMatrix = as.dist(dist_matrix), method = "complete")
sil_result_coef <- silh_result$coef
sil_result_row <- t(sil_result_coef$silIndex)

# Optimal Clusters using Silhouette Statistics
optimal_clusters_sil <- silh_result$nClus
cat("Optimal Number of Clusters by Silhouette Statistic is:", optimal_clusters_sil)

plot(1: max_clusters, sil_result_row, type = "b",
     xlab = "No of Clusters", ylab = "Silhouette Scores",
     main = "UME Silhouette Statistic Clustering Analysis")


# Gap Statistic Calculation
sc1 <- hclust(as.dist(dist_matrix), method = "complete")
memb <- cutree(sc1, 1: max_clusters)
gap_result <- gap.clus(DistanceMatrix = dist_matrix, Clusters = memb, B=100)
gap_result_ngroups <- gap_result$gap.values
gap_result_row <- t(gap_result_ngroups$gap)

# Optimal Clusters using Gap Statistics
optimal_clusters_gap <- gap_result$optim.k
cat("Optimal Number of Clusters by Gap Statistic is:", optimal_clusters_gap)

plot(1: (max_clusters-1), gap_result_row, type = "b",
     xlab = "No of Clusters", ylab = "Gap Values",
     main = "UME Gap Statistic Clustering Analysis")

# Plotting the dendrogram
par(mar = c(5, 4, 4, 2) + 0.1)
plot(sc1, labels = colnames(growth_rates), sub = "", xlab = "", cex = 0.5, las = 2)
mtext("as.dist(Macf)", side = 1, line = 1, cex = 0.8)  # First line of text
mtext('hclust(*, "complete")', side = 1, line = 2, cex = 0.8)  # Second line of text
rect.hclust(sc1, k = 3, border = "purple")

# Define years corresponding to your time points
years <- seq(2000, by = 0.25, length.out = ncol(growth_rates))

# Plot each cluster
par(mfrow = c(1, 3), mar = c(4, 4, 2, 1))
memb2 <- cutree(sc1, k=3)
for (cluster_num in 1:3) {
  cluster_indices <- which(memb2 == cluster_num)
  cluster_series <- growth_rates[cluster_indices, ]
  
  # Plot the time series for the current cluster
  matplot(years, t(cluster_series), type = "l", lty = 1, lwd = 1,col = 1:length(cluster_indices),
          main = paste("Cluster", cluster_num), xlab = "Year", ylab = "Rate")
}

cluster_counts <- table(memb2)

# Display the counts
print(cluster_counts)

# If you want to loop through and print more explicitly
for (cluster_num in names(cluster_counts)) {
  cat("Cluster", cluster_num, "has", cluster_counts[cluster_num], "series.\n")
}

for (cluster_num in unique(memb2)) {
  series_in_cluster <- names(memb2[memb2 == cluster_num])  # Get series names for this cluster
  cat("Cluster", cluster_num, "contains the following series:\n")
  print(series_in_cluster)
}
# legend("topright", legend = paste("Series", cluster_1_indices), col = 1:ncol(cluster_1_series), lty = 1)

# Jump Method Implementation
J_0 <- J_original_distribution(growth_rates)

J_k <- J_reference_bootstrap_distribution(growth_rates, B = 100)

results <- Jump_test(J_0, J_k)

if (results$conclusion == "Reject H0: There is evidence of more than one cluster (Ha: G > 1)"){
  optimal_no_of_clust <- determine_optimal_no_of_clusters(growth_rates, Gmax=max_clusters)
}else{
  optimal_no_of_clust <- 1
}

cat("Test Statistic Values:", results$T_values, "\n")
cat("Critical Values:", results$critical_values, "\n")
cat("Conclusion:", results$conclusion, "\n")
cat("Optimal Number of Clusters:", optimal_no_of_clust, "\n")
