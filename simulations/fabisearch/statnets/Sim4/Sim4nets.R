### Simulation parameters ###
reps = 1:100
simnum = as.numeric(substr(getwd(), nchar(getwd()), nchar(getwd())))
#############################

# Load all shared simulation functions
filesrcs = list.files(path = "../../../sharedfunctions", pattern="*.R")
sapply(paste0("../../../sharedfunctions/", filesrcs), source, .GlobalEnv)
rm(filesrcs)

# Run the simulations
source("../../scripts/SimNets.R")