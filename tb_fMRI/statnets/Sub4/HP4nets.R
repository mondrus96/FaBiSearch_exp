# Define the subject number based on parent folder
subnum = as.numeric(substr(getwd(), nchar(getwd()), nchar(getwd())))

# Run the network estimation script
source("../../scripts/HPnets.R")