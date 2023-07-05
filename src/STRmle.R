# Title        : Template script for the MLE of haplotype frequencies, MOI, and prevalence
#                from example dataset.
# Objective    : Estimate haplotype frequencies, MOI, and prevalence
#                from genomic/molecular data (microsatellite data)
# Created by   : Christian Tsoungui Obama
# Created on   : 05.05.22
# Last modified: 29.06.23

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
mle(data, c(2,3))

# Finding Haplotype frequencies assuming MOI known, i.e., lambda=0.2
mle(data, c(2,3), plugin = 0.2)

# Finding MLEs (haplotype frequencies and MOI) with bootstrap bias-correction ('Bootstrap')
mle(data, c(2,3), plugin = NULL, BC = TRUE, method = "bootstrap")

# Finding MLEs (haplotype frequencies and MOI) with bootstrap bias-correction ('jackknife') with plugin
mle(data, c(2,3), plugin = 0.2, BC = TRUE, method = "jackknife")

# Finding MLEs (haplotype frequencies and MOI) using a 95% confidence interval
mle(data, c(2,3), plugin = NULL, CI = TRUE)

# Finding MLEs (haplotype frequencies and MOI) using a 90% confidence interval and 20000 bootstrap samples
mle(data, c(2,3), plugin = NULL, CI = TRUE, B = 20000, alpha = 0.1)

# Finding LD between the two loci. The function outputs the LD measures D', r-squared, Q-star, and ALD.
ld(data, c(2,3), CI = TRUE, B = 15000, alpha = 0.1)

# Finding MLEs from unformatted dataset (haplotype frequencies and MOI)
########################################################################
## Step 1: Transforming the dataset using the data_format() function
##         The transformed dataset is accessed in data[[1]] 
##         and the genetic architecture is accessed in data[[3]]
########################################################################
data <- data_format(data_unf, output.id = TRUE)

########################################################################
## Step 2: Finding the MLEs for any pair of marker of interest
########################################################################

### Selecting Markers of interest (if more than two markers). Include 1 in the vector
### list to include the ID column, e.g., c(2,3) without ID, and c(1, 2, 3) with ID.
markers <- c(3,4)

### Estimating the MLEs
mle(data[[1]][,(markers+1)], data[[3]][markers], id = FALSE)

# Finding LD from unformatted dataset (haplotype frequencies and MOI)
########################################################################
## Transform the dataset using the data_format() function, 
## and drop missing information like above,
## then find LD
########################################################################
ld(data[[1]][,markers], data[[3]][markers], id = FALSE, CI = TRUE, B = 15000, alpha = 0.1)
