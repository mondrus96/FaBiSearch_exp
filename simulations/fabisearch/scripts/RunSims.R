# This function takes as an input which iterations and which simulation type and 
# returns results from running the fabisearch change point detection method

# Load library
library(fabisearch)

# Loop through and generate results
all.results = c()
for(i in reps){
  # Print the iteration
  print(i)
  
  # Set the seed
  set.seed(i*123)
  
  # Generate simulation depending on the sim type
  y = do.call(paste0("sim", simnum), list(i))
  
  # Run fabisearch
  result = detect.cps(y + 100, alpha = 0.05, testtype = "t-test", mindist = 50, nruns = 100, nreps = 1000)
  
  # Append to the overall results
  result.df = cbind(i, result$rank, result$change_points, as.numeric(result$compute_time))
  colnames(result.df)[c(2, 5)] = c("rank", "mins")
  all.results = rbind(all.results, result.df)
  
  # Print and save results
  print(all.results)
  save(all.results, file = paste0("sim", simnum, "_", min(reps), "to", max(reps), ".rda"))
}