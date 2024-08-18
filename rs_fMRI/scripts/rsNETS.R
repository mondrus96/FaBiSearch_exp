# This function takes as an input which scanning session and
# returns results from running the fabisearch network estimation

# Load library
library(fabisearch)

# Load data
data = read.table(paste0("../../data/gordon_visit_", visit, ".txt"))

# Load change point data
load(paste0("../../cpd/Visit", visit, "/gord", visit, ".rda"))
cpt.data = all.results
rm(all.results)

# Set seed and loop through subjects
set.seed(12345)
all.results = list()
for(i in 1:25){
  # Define the indices for the current subject
  indices = (((i-1)*197)+1):(i*197)
  curr.subj = data[indices, 2:ncol(data)]
  
  # Define the current subject's changepoints and rank
  curr.cpt.data = cpt.data[cpt.data$stat_test & cpt.data$i == i,]
  changepoints = curr.cpt.data$T
  rank = curr.cpt.data$rank[1]
  
  # Run fabisearch
  result = est.net(curr.subj + 100, lambda = 7, nruns = 100,
                   changepoints = changepoints, rank = rank)
  
  # Append to the overall results
  all.results[[i]] = result
  
  # Save results
  save(all.results, file = paste0("gord", visit, "_nets.rda"))
}
