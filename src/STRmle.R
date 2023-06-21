# Title        : Template script for the MLE of haplotype frequencies, MOI, and prevalence
#                from example dataset.
# Objective    : Estimate haplotype frequencies, MOI, and prevalence
#                from genomic/molecular data (microsatellite data)
# Created by   : Christian Tsoungui Obama
# Created on   : 05.05.22
# Last modified: 14.03.23

# Install the necessary packages if necessary
install.packages('openxlsx')   # Comment this line if openxlsx installed

# Loading libraries
library(openxlsx)

path <- "/Users/christian/Library/CloudStorage/GoogleDrive-christian.tsoungui@aims-cameroon.org/.shortcut-targets-by-id/1Ulru-DjbFRaMVB7Vj9tJ4NfyzPDkhzOr/Maths against Malaria/Christian/Models/MultiAllelicBiLociModel"

# Importing the dataset (formatted)
data <- read.xlsx(paste0(path,'/datasets/example_dataset1.xlsx'), 1)

# Importing the dataset (unformatted)
data_unf <- read.xlsx(paste0(path,'/datasets/example_dataset4.xlsx'), 1)

# Load external resources
source(paste0(path,'/src/STRModel.R'))#("/home/janedoe/Documents/src/STRmodel.R")

# Finding MLEs from formatted dataset (haplotype frequencies and MOI)
mle(data, c(3,2), id = TRUE, plugin = NULL)

# Finding MLEs from unformatted dataset (haplotype frequencies and MOI)
## Step 1: Transforming the dataset using the data_format() function
data <- data_format(data_unf, id=TRUE)

## Step 2: Finding the MLEs for any pair of marker of interest
### Dropping Missing data and counting alleles from 0
pick <- rowSums(data[[1]] == 0) > 0 
data[[1]] <- data[[1]][!pick,] - 1
markers <- c(1,2) # pair of markers of interest
### Estimating the MLEs
mle(data[[1]][,markers], data[[3]][markers], id = FALSE, plugin = NULL)

# Finding Haplotype frequencies assuming MOI known, i.e., lambda=0.2
mle(data, c(3,2), id = TRUE, plugin = 0.2)

# Finding MLEs (haplotype frequencies and MOI) with bootstrap bias-correction
mle(data, c(3,2), id = TRUE, plugin = NULL, BC = TRUE, method = "bootstrap")

# Finding MLEs (haplotype frequencies and MOI) using a 95% confidence interval
mle(data, c(3,2), id = TRUE, plugin = NULL, CI = TRUE)

# Finding MLEs (haplotype frequencies and MOI) using a 90% confidence interval and 15000 bootstrap samples
mle(data, c(3,2), id = TRUE, plugin = NULL, CI = TRUE, B = 15000, alpha = 0.1)

# Finding LD between the two loci. The function outputs the LD measures D', r-squared, Q-star, and ALD.
ld(data, c(3,2), id = TRUE, CI = TRUE, B = 15000, alpha = 0.1)
