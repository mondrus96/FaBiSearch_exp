# This function takes as an input which iterations and which simulation type and 
# returns results from running the fabisearch change point detection method

# Load library
library(fabisearch)

# Load data
data = read.table(paste0("../../data/HP", subnum, ".txt"))
data = data[, 2:ncol(data)]

# Load library
library(fabisearch)

# Set seed and loop through subjects
set.seed(12345)
all.results = c()
for(i in 1:4){
  # Define the current run
  curr.run = data[data$run == i, 2:ncol(data)]
  
  # Run fabisearch
  result = detect.cps(curr.run[,2:ncol(curr.run)] + 100, alpha = 0.05, testtype = "t-test", mindist = 50, nruns = 100, nreps = 1000)
  
  # Adjust change points to true time points
  result$change_points$T = curr.run$Time[result$change_points$T]
  
  # Append to the overall results
  result.df = cbind(i, result$rank, result$change_points)
  colnames(result.df)[2] = "rank"
  all.results = rbind(all.results, result.df)
  
  # Print and save results
  print(all.results)
  save(all.results, file = paste("HPcpd", subnum, ".rda", sep = ""))
}