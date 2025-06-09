# Required Libraries
library(forecast)
library(cluster)
library(stats)
library(SLBDD)
library(zoo)
library(readxl)
library(tidyr)
library(dplyr)
library(MASS)
source("JumpMethodTaiwan.R")

# Determine where data starts (e.g., after "@data")
lines <- readLines("C:\\Users\\devas\\Desktop\\MA-900\\Extra Datasets\\fred_md_dataset.tsf")
# Find the starting point of the data
data_start <- grep("@data", lines) + 1

# Extract only the time series lines
time_series_lines <- lines[data_start:length(lines)]

# Parse the time series data
parsed_data <- lapply(time_series_lines, function(line) {
  # Split by colons
  parts <- strsplit(line, ":")[[1]]
  
  # Extract the series ID and factor (e.g., temperature)
  series_id <- parts[1]

  # Extract the values (from the last part, split by commas)
  values <- as.numeric(unlist(strsplit(parts[length(parts)], ",")))
  
  # Generate timestamps for the dataset
  start_time <- as.POSIXct("1959-01-01 00:00:00", tz = "UTC")
  timestamps <- seq(from = start_time, by = "month", length.out = length(values))
  
  # Return as a data frame
  data.frame(
    time = timestamps,  # Time stamps (1, 2, 3, ...)
    series_id = series_id,
    value = values,
    stringsAsFactors = FALSE
  )
})

# Combine all parsed data into a single data frame
long_data <- do.call(rbind, parsed_data)

# Reshape data to wide format: time as rows, factors as columns
wide_data <- long_data %>%
  pivot_wider(names_from = series_id, values_from = value)

# Checking for NA values
anyNA(wide_data)

is.data.frame(wide_data)
fred_data <- as.matrix(wide_data[, -1])

# Plot the Data
custom_colors <- c("blue", "red", "purple", "orange", "brown", "pink", "black", "yellow", "cyan")
col_vector <- rep(custom_colors, length.out = ncol(fred_data))

# Generate a sequence of hourly timestamps
start_time <- as.POSIXct("1959-01-01 00:00:00", tz = "UTC")
timestamps <- seq(from = start_time, by = "month", length.out = nrow(fred_data))

par(mfrow = c(1, 1))  # Reset to single plotting
matplot(x = timestamps,
        y = fred_data, 
        type = "l",
        lty = 1,                      # Line type
        col = col_vector,  # Assign colors to series
        xlab = "Time (Years)",
        ylab = "Economic Indicators",
        main = "Federal Reserve Bank Macro-Economic Indicators"
)

fred_data[fred_data <= 0]

# Shift the dataset column-wise
shifted_data <- data.frame(matrix(ncol = ncol(fred_data), nrow = nrow(fred_data)))  # Initialize an empty data frame
colnames(shifted_data) <- colnames(fred_data)  # Preserve column names

fred_data <- as.data.frame(fred_data)
for (col_name in colnames(fred_data)) {
  col <- (fred_data[[col_name]])  # Extract the column
  shift_value <- abs(min(col, na.rm = TRUE)) + 1e-6  # Compute shift value
  shifted_data[[col_name]] <- col + shift_value  # Apply shift and store in the new dataset
}

shifted_data[shifted_data <=0]

# Plotting Shifted FRED Data
par(mfrow = c(1, 1))  # Reset to single plotting
matplot(x = timestamps,
        y = shifted_data, 
        type = "l",
        lty = 1,                      # Line type
        col = col_vector,  # Assign colors to series
        xlab = "Time (Years)",
        ylab = "Economic Indicators",
        main = "Shifted Federal Reserve Bank Macro-Economic Indicators"
)

# Taking Log Returns of the Shifted Data
shifted_data <- as.matrix(shifted_data)
log_returns_FRED <- diff(log(shifted_data))

# Plotting Log Returns of FRED Data
timestamps <- timestamps[-1]
par(mfrow = c(1, 1))  # Reset to single plotting
matplot(x = timestamps,
        y = log_returns_FRED, 
        type = "l",
        lty = 1,                      # Line type
        col = col_vector,  # Assign colors to series
        xlab = "Time (Years)",
        ylab = "Economic Indicators",
        main = "Log Returns of Federal Reserve Bank Macro-Economic Indicators"
)

num_series <- ncol(log_returns_FRED)

acf_features <- matrix(NA, nrow = num_series, ncol = 6)
for (i in 1: num_series) {
  acf_result <- acf(log_returns_FRED[, i], plot = FALSE, lag.max = 6)
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
plot(sc1, labels = colnames(log_returns_FRED), sub = "", xlab = "", main = "logDiffFredMD ACF ",cex = 0.5, las = 2)
rect.hclust(sc1, k = 4, border = "blue")
mtext("dist((MACF))\n hclust(*, 'complete')", side = 1, line = 4, cex = 0.8)


# Plot each cluster
timestamps <- seq(from = start_time, by = "month", length.out = ncol(log_returns_FRED))
par(mfrow = c(1, 3), mar = c(4, 4, 2, 1))
memb2 <- cutree(sc1, k=as.numeric(optimal_clusters_gap))
for (cluster_num in 1:as.numeric(optimal_clusters_gap)) {
  cluster_indices <- which(memb2 == cluster_num)
  cluster_series <- log_returns_FRED[cluster_indices, ]
  
  # Plot the time series for the current cluster
  matplot(timestamps, t(cluster_series), type = "l", lty = 1, lwd = 1, col = 1:length(cluster_indices),
          main = paste("Cluster", cluster_num), xlab = "Year", ylab = "Economic Indicators Log Returns")
}

# Jump Method Implementation
optimal_no_of_clust <- 0
J_0 <- J_original_distribution(log_returns_FRED)

J_k <- J_reference_bootstrap_distribution(log_returns_FRED, B = 100)

results <- Jump_test(J_0, J_k)

if (results$conclusion == "Reject H0: There is evidence of more than one cluster (Ha: G > 1)"){
  optimal_no_of_clust <- determine_optimal_no_of_clusters(log_returns_FRED, Gmax=max_clusters)
}else{
  optimal_no_of_clust <- 1
}

cat("Test Statistic Values:", results$T_values, "\n")
cat("Critical Values:", results$critical_values, "\n")
cat("Conclusion:", results$conclusion, "\n")
cat("Optimal Number of Clusters:", optimal_no_of_clust, "\n")
