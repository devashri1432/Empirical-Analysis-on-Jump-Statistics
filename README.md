# Empirical-Analysis-on-Jump-Statistics
Conducted research on clustering univariate scalar time series using the Jump Statistic, a novel method based on dendrogram jumps and AR-sieve bootstrapping. Compared performance with Silhouette and Gap Statistics using Monte Carlo simulations and real-world datasets (EDEN ISS, FRED-MD).

This research study focuses on clustering univariate scalar time series to uncover hidden patterns and identify groups with similar temporal properties. The project addresses two fundamental challenges in time series clustering:

**_1. Determining the presence of multiple clusters
2. Accurately estimating the number of clusters**
_
Three statistical methods—_Jump Statistic, Silhouette Statistic_, and _Gap Statistic_—are applied in conjunction with agglomerative hierarchical clustering. The primary contribution is a detailed empirical evaluation of the Jump Statistic, a relatively new technique with limited real-world use. Key innovations include:

**1. Designing a bootstrap-based reference distribution using the AR-sieve method
2. Performing hypothesis testing via upper quantile comparison of dendrogram jumps
3. Implementing a novel iterative algorithm to determine the optimal cluster count
4. Validating the method via Monte Carlo simulations and empirical datasets (EDEN ISS and FRED-MD)**

This thesis demonstrates the effectiveness of the Jump method in high-dimensional time series scenarios, outperforming traditional metrics in various simulated and real-world cases.
