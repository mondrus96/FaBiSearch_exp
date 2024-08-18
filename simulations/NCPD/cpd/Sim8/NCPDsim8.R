# Define the cluster number to use
try.K = 2

# Get simulation number
simnum = as.numeric(substr(getwd(), nchar(getwd()), nchar(getwd())))

# Load all shared simulation functions
filesrcs = list.files(path = "../../../sharedfunctions", pattern="*.R")
sapply(paste0("../../../sharedfunctions/", filesrcs), source, .GlobalEnv)
rm(filesrcs)

# Load all the functions for NCPD
source("../../functions/SpectralFunctions.R")
source("../../functions/ParallelFunctions.R")

# Run the experiment
source("../../functions/RunSimsNCPD.R")