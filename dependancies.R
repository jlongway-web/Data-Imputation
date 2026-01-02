# dependencies.R
# Run this script to install all necessary libraries for the project

install.packages(c(
  "car",      # For regression diagnostics
  "ggplot2",  # For high-quality visualizations
  "GGally",   # For correlation matrices
  "dplyr",    # For data manipulation
  "tidyr",    # For data reshaping (pivoting)
  "modeest",  # For calculating the mode
  "olsrr",    # For stepwise regression
  "kernlab",  # For SVM models
  "pROC"      # For ROC and AUC analysis
))

print("All dependencies have been installed.")
