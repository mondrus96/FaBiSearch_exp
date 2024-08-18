# This function takes as an input which iterations and which simulation type and 
# returns results from running the fabisearch change point detection method

# Load library
library(fabisearch)

# Define breaks and true.K, which contains the indices of stationary segments
# and true cluster number for those segments
if(simnum == 1){
  true.K = 2
  breaks = list(1:200)
} else if (simnum == 2){
  true.K = c(2, 2)
  breaks = list(1:100, 101:200)
} else if (simnum == 3){
  true.K = c(3, 2, 2, 3)
  breaks = list(1:100, 101:200, 201:300, 301:400)
} else if (simnum == 4){
  true.K = c(2, 2, 2)
  breaks = list(1:200, 201:400, 401:600)
} else if (simnum == 5){
  true.K = c(2, 2, 2)
  breaks = list(1:100, 101:200, 201:300)
} else if (simnum == 6){
  true.K = c(7, 2, 7)
  breaks = list(1:100, 101:200, 201:300)
} else if (simnum == 7){
  true.K = c(2, 2)
  breaks = list(1:100, 151:200)
} else if (simnum == 8){
  true.K = c(2, 2)
  breaks = list(1:100, 106:200)
}

# Load results, loop through and append the separate runs
filesrcs = list.files(path = paste0("../../cpd/Sim", simnum), pattern="*.rda")
cpt.results = c()
for(i in 1:length(filesrcs)){
  # Load the current run
  load(paste0("../../cpd/Sim", simnum, "/", filesrcs[i]))
  
  # Append the results
  cpt.results = rbind(cpt.results, all.results)
  
  # Remove the loaded run
  rm(all.results)
}
rm(filesrcs)

# Loop through
all.results = list()
for(i in reps){
  # Print the iteration
  print(i)
  
  # Generate simulation depending on the sim type
  y = do.call(paste0("sim", simnum), list(i))
  
  # Narrow down to final change points, append beginning and end
  curr.cpt = c(0, cpt.results$T[cpt.results$stat_test & cpt.results$i == i], nrow(y))
  
  # Find the rank used as optimal
  rank = cpt.results$rank[cpt.results$i == i][1]
  
  # Loop through the stationary segments
  curr.rep = list()
  for(j in 1:(length(curr.cpt) - 1)){
    # Define the current indices based on detected change point
    curr.indices = (curr.cpt[j] + 1):(curr.cpt[j + 1])
    
    # Define the current stationary segment
    curr.seg = y[curr.indices,]
    
    # See where the greatest overlap is, use this to define lambda
    overlap = c()
    for(k in 1:length(breaks)){
      overlap[k] = sum(curr.indices %in% breaks[[k]])
    }
    lambda = true.K[match(max(overlap), overlap)]
    
    # Run network estimation
    result = est.net(curr.seg + 100, lambda = lambda, nruns = 100, rank = rank)
    
    # Append to the current repetition results
    curr.rep[[j]] = result
  }
  # Append to all results
  all.results[[i]] = curr.rep
  
  # Save results
  save(all.results, file = paste0("sim", simnum, "_nets.rda"))
}
