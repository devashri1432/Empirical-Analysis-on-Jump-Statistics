# Required Libraries
library(forecast)
library(cluster)
library(stats)
library(SLBDD)
library(zoo)
library(readxl)
library(tidyr)
library(dplyr)
library(data.table)
source("JumpMethodTaiwan.R")

# Define the parent directory containing the folders
parent_dir <- "C:\\Users\\devas\\Desktop\\MA-900\\Extra Datasets\\Telemetry Data 2020\\edeniss2020"

# Get the paths of all CSV files from all folders
csv_files <- list.files(parent_dir, pattern = "\\.csv$", full.names = TRUE, recursive = TRUE)

# Initialize an empty list to store datasets
list_of_data <- list()

# Flag to track if the time column is already added
time_column_added <- FALSE
time_column <- NULL

# Process each file
for (file in csv_files) {
  data <- fread(file)  # Read the entire file
  
  # Extract the time column only once (assuming it's named "time")
  if (!time_column_added && "time" %in% names(data)) {
    time_column <- data[, .(time)]  # Retain the time column
    time_column_added <- TRUE      # Mark that the time column has been added
  }
  
  # Remove the time column from the current dataset (if it exists)
  data <- data[, !("time"), with = FALSE]
  
  # Add the remaining data to the list
  list_of_data <- append(list_of_data, list(data))
}

# Combine all datasets column-wise
feature_data <- do.call(cbind, list_of_data)

# Set the column names to avoid file name prefixes
colnames(feature_data) <- unlist(lapply(list_of_data, names))

# Combine the single time column with the feature data
sensor_data <- data.frame(time_column, feature_data)

# View the final combined dataframe
head(sensor_data)

is.data.frame(sensor_data)

sensor_data <- as.data.frame(sensor_data)

colnames(sensor_data)

cols_to_drop <- c("Filename", "Path", "Subsystem", "Sensor.Type..short.", "Sensor.Type..long.", "Unit")

sensor_data <- sensor_data[, !names(sensor_data) %in% cols_to_drop]
head(sensor_data)

required_data <- sensor_data[, -c(1)]
anyNA(required_data)

is.numeric(required_data)
numeric_columns <- sapply(required_data, is.numeric)
print(numeric_columns)

summary(required_data)

custom_colors <- c("blue", "red", "green", "purple", "orange", "brown", "pink", "cyan", "black", "yellow")
col_vector <- rep(custom_colors, length.out = ncol(required_data))

sensor_matrix <- as.matrix(required_data)

start_time <- as.POSIXct("2020-01-01 00:05:00", tz = "UTC")
timestamps <- seq(from = start_time, by = "min", length.out = nrow(required_data))

par(mfrow = c(1, 1))  # Reset to single plotting
matplot(
  x = timestamps,  # Index for rows
  y = sensor_matrix,  # Convert data to a matrix
  type = "l",  # Line plot
  col = 1:ncol(required_data),  # Unique color for each column
  lty = 1,  # Line type
  xlab = "Time (Min)",
  ylab = "Sensor Readings",
  main = ""
)

# Compute rolling dynamic quantiles
window_size <- 3  # Rolling window size (e.g., 24 hours)
rolling_quantiles <- rollapply(
  sensor_data,
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
  x = 1:nrow(sensor_data), 
  y = sensor_data, 
  type = "l", 
  lty = 1,
  col = "black",
  xlab = "Time (Hours)", 
  ylab = "PM-2.5 Values", 
  main = ""
)

lines(1:nrow(sensor_data), q5, col = "red", lty = 1, lwd = 1)  # 25th percentile
lines(1:nrow(sensor_data), q50, col = "blue", lty = 1, lwd = 1)  # Median
lines(1:nrow(sensor_data), q95, col = "green", lty = 1, lwd = 1)  # 75th percentile
# par(mfrow = c(5, 5), mar = c(2, 2, 2, 2))  # Reduce margins
# 
# # Generate ACF plots for each time series
# for (i in 1:ncol(required_data)) {
#   acf(required_data[[i]], main = paste("ACF of TS", i), lag.max = 50)
# }
# 
# # Reset layout
# par(mfrow = c(1, 1))
# 
# diff_sensor_data <- diff(as.matrix(required_data), lag = 288)
# par(mfrow = c(5, 5), mar = c(2, 2, 2, 2))  # Reduce margins
# 
# # Generate ACF plots for each time series
# for (i in 1:ncol(diff_sensor_data)) {
#   acf(diff_sensor_data[[i]], main = paste("ACF of TS", i), lag.max = 10)
# }
# 
# diff_sensor_data <- as.data.frame(diff_sensor_data)
# 
# acf(diff_sensor_data$par.1, lag.max = 10)
# 
# # Reset layout
# par(mfrow = c(1, 1))
# 
# decomposed <- stl(ts(diff_sensor_data[[1]], frequency = 12), s.window = "periodic")
# plot(decomposed)

nsdiffs(sensor_matrix)
ndiffs(sensor_matrix)

sensor_data <- as.data.frame(sensor_matrix)

seasonal_diffs <- sapply(sensor_data, ndiffs)
print(seasonal_diffs)

diff_sensor_data <- diff(sensor_matrix)
timestamps <- timestamps[-1]
par(mfrow = c(1, 1))  # Reset to single plotting
matplot(
  x = timestamps,  # Index for rows
  y = diff_sensor_data,  # Convert data to a matrix
  type = "l",  # Line plot
  col = 1:ncol(diff_sensor_data),  # Unique color for each column
  lty = 1,  # Line type
  xlab = "Time (Min)",
  ylab = "Differenced Sensor Readings",
  main = ""
)


par(mfrow = c(1, 1))  # Reset to single plotting
matplot(
  x = timestamps,  # Index for rows
  y = full_diff_data,  # Convert data to a matrix
  type = "l",  # Line plot
  col = 1:ncol(full_diff_data),  # Unique color for each column
  lty = 1,  # Line type
  xlab = "Time (Min)",
  ylab = "Differenced Sensor Readings",
  main = ""
)

num_series <- ncol(diff_sensor_data)

acf_features <- matrix(NA, nrow = num_series, ncol = 6)
for (i in 1: num_series) {
  acf_result <- acf(diff_sensor_data[, i], plot = FALSE, lag.max = 6)
  acf_features[i, ] <- acf_result$acf[2:7]
}

max_clusters <- 5
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
plot(sc1, labels = colnames(diff_sensor_data), sub = "", xlab = "", main = "diffEdenSensor ACF ",cex = 0.5, las = 2)
mtext("as.dist(Macf)", side = 1, line = 1, cex = 0.8)  # First line of text
mtext('hclust(*, "complete")', side = 1, line = 2, cex = 0.8)  # Second line of text
rect.hclust(sc1, k = 3, border = "blue")

# Plot each cluster
timestamps <- seq(from = start_time, by = "min", length.out = ncol(diff_sensor_data))
par(mfrow = c(1, 3), mar = c(4, 4, 2, 1))
memb2 <- cutree(sc1, k=5)
for (cluster_num in 1:5) {
  cluster_indices <- which(memb2 == cluster_num)
  cluster_series <- diff_sensor_data[cluster_indices, ]
  
  # Plot the time series for the current cluster
  matplot(timestamps, t(cluster_series), type = "l", lty = 1, lwd = 1, col = 1:length(cluster_indices),
          main = paste("Cluster", cluster_num), xlab = "Time", ylab = "EDEN ISS Sensor Readings")
}

# Jump Method Implementation
optimal_no_of_clust <- 0
J_0 <- J_original_distribution(diff_sensor_data)

J_k <- J_reference_bootstrap_distribution(diff_sensor_data, B = 100)

results <- Jump_test(J_0, J_k)

if (results$conclusion == "Reject H0: There is evidence of more than one cluster (Ha: G > 1)"){
  optimal_no_of_clust <- determine_optimal_no_of_clusters(diff_sensor_data, Gmax=max_clusters)
}else{
  optimal_no_of_clust <- 1
}

cat("Test Statistic Values:", results$T_values, "\n")
cat("Critical Values:", results$critical_values, "\n")
cat("Conclusion:", results$conclusion, "\n")
cat("Optimal Number of Clusters:", optimal_no_of_clust, "\n")

# Get cluster assignments
memb2 <- cutree(sc1, k = 5)

# Ensure clusters are correctly ordered
unique_clusters <- sort(unique(memb2))  

# Initialize matrix for storing average ACF values
num_clusters <- length(unique_clusters)
avg_acf_clusters <- matrix(NA, nrow = num_clusters, ncol = 6)  # 7 clusters, 6 lags
cluster_sizes <- numeric(num_clusters)  # Store cluster sizes

# par(mfrow = c(1, 3), mar = c(4, 4, 2, 1))
# memb2 <- cutree(sc1, k=5)
# for (cluster_num in 1:5) {
#   cluster_indices <- which(memb2 == cluster_num)
#   cluster_series <- required_data[cluster_indices, ]
#   
#   # Plot the time series for the current cluster
#   matplot(years, t(cluster_series), type = "l", lty = 1, lwd = 1,col = 1:length(cluster_indices),
#           main = paste("Cluster", cluster_num), xlab = "Year", ylab = "Rate")
# }

# Loop over each cluster
for (cluster_num in unique_clusters) {
  # Get indices of time series (columns) belonging to the current cluster
  cluster_indices <- which(memb2 == cluster_num)
  
  # Extract time series for this cluster (Columns: Time Series, Rows: Time Points)
  cluster_series <- diff_sensor_data[, cluster_indices, drop = FALSE]
  
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

colnames(diff_sensor_data)[68]

colnames(diff_sensor_data)[c(67,69,70,72)]

colnames(diff_sensor_data)[c(1,9,10,11,50,53,75,77,84,97)]

colnames(diff_sensor_data)[c(2,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,
                             39,40,41,42,43,44,45,46,47,48,49,51,52,54,55,88,90,91,92,93,94)]


print(names(memb2))

ncol(sensor_data)
decomposed_ts <- decompose(ts(sensor_data$co2.1, frequency = 7), type="additive")
plot(decomposed_ts)


# Plot individual time series
par(mfrow = c(3,3))

for (i in 1:ncol(required_data)){
  plot(required_data[,i], type="l",
       main = colnames(required_data)[i],
       xlab = "Time", ylab = "Value")
}

plot(required_data$temp.iocs_out)
plot(decompose(ts(required_data$temp.iocs_out, frequency = 12))
)
nsdiffs(required_data$temp.iocs_out)
ndiffs(required_data$r3.4l)

nsdiffs(as.data.frame(diff_sensor_data)$temp.iocs_out)

par(mfrow = c(3,4))
for (i in 1:ncol(required_data)){

decomp <- decompose(ts(required_data[,i], frequency = 288))
zoom_range <- 1:1000
plot(zoom_range, decomp$trend[zoom_range], type = "l", col = "blue",
     main = colnames(required_data)[i], xlab = "Time", ylab = "Seasonality")
}

"Seasonal Differencing"

seasonal_diff_data <- apply(required_data,2, function(x) diff(x, lag = 288))
full_diff_data <- apply(seasonal_diff_data, 2, diff)

adf_results <- apply(full_diff_data, 2, function(x) {
  # ADF test can throw errors if there's not enough data, so use tryCatch
  res <- tryCatch({
    test <- adf.test(x)
    c(p.value = test$p.value, statistic = test$statistic)
  }, error = function(e) {
    c(p.value = NA, statistic = NA)
  })
  return(res)
})

# Format results as data frame
adf_df <- as.data.frame(t(adf_results))
adf_df$Series <- colnames(full_diff_data)
rownames(adf_df) <- NULL

# View the results
print(adf_df)

