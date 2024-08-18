# Define the subject number based on parent folder
subnum = as.numeric(substr(getwd(), nchar(getwd()), nchar(getwd())))

# Run change point detection
source("../../scripts/HPcpd.R")