# Required Libraries
library(forecast)
library(cluster)
library(stats)
library(SLBDD)
library(zoo)
source("JumpMethodTaiwan.R")

# data("TaiwanAirBox032017")
data_matrix <- as.matrix(TaiwanAirBox032017)
colnames(data_matrix) <- paste("Series", 1:ncol(data_matrix), sep = "_")

col_names <- colnames(data_matrix)

# Removing the Outliers
outliers <- c(1, 29, 35, 46, 70, 118, 155, 157)

remain_data <- data_matrix[, -outliers]
ncol(remain_data)

custom_colors <- c("blue", "red", "green", "purple", "orange", "brown", "pink", "cyan", "black", "yellow")
col_vector <- rep(custom_colors, length.out = ncol(remain_data))

start_date <- as.POSIXct("2017-03-01", format = "%Y-%m-%d")  # Convert to POSIXct
hours <- seq(from = start_date, by = "hour", length.out = nrow(remain_data))

matplot(x = hours,
        y = remain_data, 
        type = "l",
        lty = 1,                      # Line type
        col = col_vector,  # Assign colors to series
        xlab = "Time (Hours)",
        ylab = "PM-2.5 Values",
        main = ""
)

# Compute rolling dynamic quantiles
window_size <- 3  # Rolling window size (e.g., 24 hours)
rolling_quantiles <- rollapply(
  remain_data,
  width = window_size,
  FUN = function(x) quantile(x, probs = c(0.05, 0.5, 0.95), na.rm = TRUE),
  by.column = FALSE,
  fill = NA,
  align = "center"
)

# Separate quantiles
q5 <- rolling_quantiles[, 1]  # 25th percentile
q50 <- rolling_quantiles[, 2]  # Median
q95 <- rolling_quantiles[, 3]  # 75th percentile

matplot(
  x = 1:nrow(remain_data), 
  y = remain_data, 
  type = "l", 
  lty = 1,
  col = "black",
  xlab = "Time (Hours)", 
  ylab = "PM-2.5 Values", 
  main = ""
)

lines(1:nrow(remain_data), q5, col = "red", lty = 1, lwd = 1)  # 25th percentile
lines(1:nrow(remain_data), q50, col = "blue", lty = 1, lwd = 1)  # Median
lines(1:nrow(remain_data), q95, col = "green", lty = 1, lwd = 1)  # 75th percentile


first_diff_data <- diff(remain_data)
num_series <- ncol(first_diff_data)

acf_features <- matrix(NA, nrow = num_series, ncol = 6)
for (i in 1: num_series) {
  acf_result <- acf(first_diff_data[, i], plot = FALSE, lag.max = 6)
  acf_features[i, ] <- acf_result$acf[2:7]
}

max_clusters <- 10
# Computing the similarity matrix by calculating the euclidean distance between each series
dist_matrix <- as.matrix(dist(acf_features, method = "euclidean"))

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
plot(sc1, labels = colnames(first_diff_data), sub = "", xlab = "", main = "diffTaiwanAirBox ACF ",cex = 0.5, las = 2)
mtext("as.dist(Macf)", side = 1, line = 1, cex = 0.8)  # First line of text
mtext('hclust(*, "complete")', side = 1, line = 2, cex = 0.8)  # Second line of text
rect.hclust(sc1, k = 4, border = "blue")

# Plot each cluster
hours <- hours[-1]
par(mfrow = c(1, 3), mar = c(4, 4, 2, 1))

# Get cluster assignments
memb2 <- cutree(sc1, k = 7)

# Ensure clusters are correctly ordered
unique_clusters <- sort(unique(memb2))  

# Initialize matrix for storing average ACF values
num_clusters <- length(unique_clusters)
avg_acf_clusters <- matrix(NA, nrow = num_clusters, ncol = 6)  # 7 clusters, 6 lags
cluster_sizes <- numeric(num_clusters)  # Store cluster sizes

# Loop over each cluster
for (cluster_num in unique_clusters) {
  # Get indices of time series (columns) belonging to the current cluster
  cluster_indices <- which(memb2 == cluster_num)
  
  # Extract time series for this cluster (Columns: Time Series, Rows: Time Points)
  cluster_series <- first_diff_data[, cluster_indices, drop = FALSE]
  
  # Compute ACF for each time series (column-wise processing)
  acf_values <- do.call(cbind, lapply(1:ncol(cluster_series), function(i) {
    acf_res <- acf(cluster_series[, i], lag.max = 6, plot = FALSE)  # Compute ACF
    acf_lags <- acf_res$acf[2:7]  # Extract first 6 lags (excluding lag 0)
    
    # Ensure it has exactly 6 values, filling with NA if needed
    if (length(acf_lags) < 6) {
      acf_lags <- c(acf_lags, rep(NA, 6 - length(acf_lags)))
    }
    
    return(acf_lags)
  }))
  
  # Compute the average ACF for this cluster
  avg_acf_clusters[which(unique_clusters == cluster_num), ] <- rowMeans(acf_values, na.rm = TRUE)
  
  # Store cluster size
  cluster_sizes[which(unique_clusters == cluster_num)] <- length(cluster_indices)
}

# Convert results to a DataFrame
acf_df <- as.data.frame(t(avg_acf_clusters))  # Transpose to match row-wise format
colnames(acf_df) <- paste0("Cluster ", unique_clusters)  # Rename columns
acf_df$Lag <- paste0("\u0304r", 1:nrow(acf_df))  # Assign labels for ACF lags

# Append cluster sizes
sizes_df <- as.data.frame(t(cluster_sizes))
colnames(sizes_df) <- colnames(acf_df)[-ncol(acf_df)]
sizes_df$Lag <- "Size"

# Combine results
final_table <- rbind(acf_df, sizes_df)

# Print table
print(final_table)



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

# Jump Method Implementation
optimal_no_of_clust <- 0
J_0 <- J_original_distribution(first_diff_data)

J_k <- J_reference_bootstrap_distribution(first_diff_data, B = 100)

results <- Jump_test(J_0, J_k)

if (results$conclusion == "Reject H0: There is evidence of more than one cluster (Ha: G > 1)"){
  optimal_no_of_clust <- determine_optimal_no_of_clusters(first_diff_data, Gmax=max_clusters)
}else{
  optimal_no_of_clust <- 1
}

cat("Test Statistic Values:", results$T_values, "\n")
cat("Critical Values:", results$critical_values, "\n")
cat("Conclusion:", results$conclusion, "\n")
cat("Optimal Number of Clusters:", optimal_no_of_clust, "\n")

