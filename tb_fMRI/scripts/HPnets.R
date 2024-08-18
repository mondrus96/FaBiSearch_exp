# This function takes as an input which iterations and which simulation type and 
# returns results from running the fabisearch change point detection method

# Load library
library(fabisearch)

# Load data
data = read.table(paste0("../../data/HP", subnum, ".txt"))
data = data[, 2:ncol(data)]

# Load change point data
load(paste0("../../cpd/Sub", subnum, "/HPcpd", subnum, ".rda"))
cpt.data = all.results
rm(all.results)

# Set seed and loop through runs
set.seed(12345)
all.results = list()
for(i in 1:4){
  # Define the current run, and keep the timepoints
  curr.run = data[data$run == i, 2:ncol(data)]
  
  # Define the current run's changepoints and rank
  curr.cpt.data = cpt.data[cpt.data$stat_test & cpt.data$i == i,]
  changepoints = which(curr.run$Time %in% curr.cpt.data$T)
  rank = curr.cpt.data$rank[1]
  
  # Run fabisearch
  result = est.net(curr.run[, 2:ncol(curr.run)] + 100, lambda = seq(0.01, 0.99, 0.01), nruns = 100,
                   changepoints = changepoints, rank = rank)
  
  # Append to the overall results
  all.results[[i]] = result
  
  # Print and save results
  save(all.results, file = paste0("HPnets", subnum, ".rda"))
}