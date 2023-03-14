# Title        : Template script for the MLE of haplotype frequencies, MOI, and prevalence
#                from example dataset.
# Objective    : Estimate haplotype frequencies, MOI, and prevalence
#                from genomic/molecular data (microsatellite data)
# Created by   : Christian Tsoungui Obama
# Created on   : 05.05.22
# Last modified: 09.03.23

# Install the necessary packages if necessary
install.packages('openxlsx')   # Comment this line if xlsx installed

# Loading libraries
library(openxlsx)

# Import the dataset
DATA <- read.xlsx('/home/janedoe/Documents/example_data1.xlsx', 1)

# Load external resources
source("/home/janedoe/Documents/STRModel.R")

# Find MLEs (haplotype frequencies and MOI)
mle(X, c(3,2), id = TRUE, plugin = NULL)

# Find Haplotype frequencies assuming MOI known, i.e., lambda=0.2
mle(X, c(3,2), id = TRUE, plugin = 0.2)

# Find MLEs (haplotype frequencies and MOI) with bootstrap bias-correction
mle(X, c(3,2), id = TRUE, plugin = NULL, BC = TRUE, method = "bootstrap")

# Find MLEs (haplotype frequencies and MOI) using a 95% confidence interval
mle(X, c(3,2), id = TRUE, plugin = NULL, CI = TRUE)

# Find MLEs (haplotype frequencies and MOI) using a 90% confidence interval and 15000 bootstrap samples
mle(X, c(3,2), id = TRUE, plugin = NULL, CI = TRUE, B = 15000, alpha = 0.1)
