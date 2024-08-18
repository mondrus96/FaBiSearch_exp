# Define the subject number based on parent folder
visit = as.numeric(substr(getwd(), nchar(getwd()), nchar(getwd())))

# Run the network estimation script
source("../../scripts/rsNETS.R")