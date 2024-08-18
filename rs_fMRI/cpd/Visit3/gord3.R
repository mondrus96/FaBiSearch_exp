# Define the subject number based on parent folder
visit = as.numeric(substr(getwd(), nchar(getwd()), nchar(getwd())))

# Run change point detection
source("../../scripts/rsCPD.R")