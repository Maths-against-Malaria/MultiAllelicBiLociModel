# Title        : Template script for the MLE of haplotype frequencies, MOI, and prevalence
#                from example dataset.
# Objective    : Estimate haplotype frequencies, MOI, and prevalence
#                from genomic/molecular data (microsatellite data)
# Created by   : Christian Tsoungui Obama
# Created on   : 05.05.22
# Last modified: 29.06.23

path <- "/Users/christian/Library/CloudStorage/GoogleDrive-christian.tsoungui@aims-cameroon.org/.shortcut-targets-by-id/1Ulru-DjbFRaMVB7Vj9tJ4NfyzPDkhzOr/Maths against Malaria/Christian/Models/MultiAllelicBiLociModel"

# Load external resources
source(paste0(path,'/src/STRModel.R'))#("/home/janedoe/Documents/src/STRmodel.R")

# Install the necessary packages if necessary
#install.packages('openxlsx')   # Comment this line if openxlsx installed

# Loading libraries
library(openxlsx)

#################################
### Import Datasets
##################################

## Case1: Dataset in standard format
# Import the dataset
DATA <- read.xlsx(paste0(path,'/datasets/example_dataset1.xlsx'), 1)

## Case2: Dataset in natural format with multiple markers
# Import the dataset
data_unf <- read.xlsx(paste0(path,'/datasets/example_dataset4.xlsx'), 1)

# Transform the data to the standard format
Ex.data <- data_format(data_unf, output.id = TRUE)

#################################
### Estimate MLEs
##################################

# Estimate MLEs
mle(DATA, c(2,3))

# Estimate MLEs (case 2)
## Choose markers of interests 
markers <- c(3,4)
mle(Ex.data[[1]][,(markers+1)], Ex.data[[3]][markers], id = FALSE)

# Estimate MLEs with plugin estimate for MOI parameter
mle(data, c(2,3), plugin = 0.2)

# Finding MLEs (haplotype frequencies and MOI) with bootstrap bias-correction ('Bootstrap')
mle(data, c(2,3), BC = TRUE, method = "bootstrap", Bbias = 15000)

# Finding MLEs (haplotype frequencies and MOI) with bootstrap bias-correction ('jackknife') with plugin
mle(data, c(2,3), plugin = 0.2, BC = TRUE, method = "jackknife")

# Finding MLEs (haplotype frequencies and MOI) using a 95% confidence interval
mle(data, c(2,3), CI = TRUE, B = 15000)

# Finding MLEs (haplotype frequencies and MOI) using a 90% confidence interval and 20000 bootstrap samples
mle(data, c(2,3),  CI = TRUE, B = 20000, alpha = 0.1)

# Finding LD between the two loci. The function outputs the LD measures D', r-squared, Q-star, and ALD.
ld(data, c(2,3), CI = TRUE)
