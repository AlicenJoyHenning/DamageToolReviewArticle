# Install robustbase if not installed
# install.packages("robustbase")

library(robustbase)
robust_dist <- covMcd(df)$mah # Robust Mahalanobis distance

# Threshold based on Chi-square distribution (p-value < 0.01)
cutoff <- qchisq(0.99, df = ncol(df)) # df is the number of variables
outliers_robust <- robust_dist > cutoff
df$outlier_robust <- outliers_robust
