# This function takes as an input the visit number and returns 
# results from running the fabisearch change point detection method

# Load library
library(fabisearch)

# Load data
data = read.table(paste0("../../data/gordon_visit_", visit, ".txt"))

# Loop through subjects
all.results = c()
for(i in 1:25){
  # Set seed
  set.seed(i*123)
  
  # Define the indices for the current subject
  indices = (((i-1)*197)+1):(i*197)
  curr.subj = data[indices, 2:ncol(data)] + 100
  
  # Run fabisearch
  result = detect.cps(curr.subj, alpha = 0.05, testtype = "t-test", mindist = 50, nruns = 100, nreps = 1000)
  
  # Append to the overall results
  result.df = cbind(i, result$rank, result$change_points)
  colnames(result.df)[2] = "rank"
  all.results = rbind(all.results, result.df)
  
  # Print and save results
  print(all.results)
  save(all.results, file = paste0("gord", visit, ".rda"))
}