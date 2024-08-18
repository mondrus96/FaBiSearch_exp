# This file runs the analysis of simulation results
source("sim_functions.R")

##### Define a list of true change points for the simulations #####
SIMtruecpts = list(NULL,
                   100,
                   c(100, 200, 300),
                   c(200, 400),
                   c(100, 200),
                   c(100, 200),
                   100:150,
                   100:105)
SIMlengths = c(200, 200, 400, 600, 300, 300, 250, 205)
###################################################################

### Organize and save results ###
# Set working directory to the results folder
setwd("../results")
# Aggregate results from across folders
FBS = FBScpts()
NCPD = NCPDcpts()
# Calculate results, save as .csv files for adding to LaTeX file
aggFBS = computeResults(FBS, "FBS")
aggNCPD = computeResults(NCPD, "NCPD")
# Aggregate results into single csv file for plotting
plotComp(aggFBS, aggNCPD)

### Run the stationary network comparison ###
netResults(FBS = FBS, SIMlengths = SIMlengths)