###############################################################################################
### This script contains functions for doing analysis for the resting-state fMRI experiment ###
###############################################################################################

# Load libraries
library(fabisearch)
library(rgl)

##### 1) Load all the data and organize final networks #####
rsnets = function(visits){
  # visits = which visits to estimate networks for
  
  # Define the final nets variable which will store the finalized networks
  allnets = list()
  # Loop through the defined visits
  for(i in visits){
    # Load results, name based on visit number
    load(paste0("../../statnets/Visit", i, "/gord", i, "_nets", ".rda"))
    
    # Append to the allnets object
    allnets[[which(i == visits)]] = all.results
  }
  # Name the components of the list
  names(allnets) = paste0("Visit", visits)
  
  # Return the allnets object
  return(allnets)
}

##### 2) Load all the change points into a data frame #####
rscpts = function(visits){
  # visits = which visits to estimate networks for
  
  # Define the output
  output = c()
  # Loop through the defined visits
  for(i in visits){
    # Load results, name based on visit number
    load(paste0("../../cpd/Visit", i, "/gord", i, ".rda"))
    
    # Save only the final changepoints
    all.results = all.results[all.results$stat_test, c(1, 3)]
    
    # Save as csv file based on visit number
    write.csv(all.results, file = paste0("gordon", i, "table.csv"), row.names = FALSE)
    
    # Save it in the output dataframe
    output = rbind(output, cbind(i, all.results))
  }
  # Change the column names
  colnames(output) = c("visit", "subject", "time")
  rownames(output) = NULL
  
  # Return the results
  return(output)
}

##### 3) Plot and save all networks #####
plotall = function(allnets, view, comms){
  # allnets = all stationary networks
  # view = view from which to plot the rgl window from
  # comms = communities to plot
  
  # Loop through visits
  for(i in 1:length(allnets)){
    print(names(allnets)[i])
    # Loop through subjects
    for(j in 1:length(allnets[[i]])){
      print(paste0("Subject ", j))
      # Loop through the networks
      for(k in 1:length(allnets[[i]][[j]])){
        print(paste0("Network ", k))
        # Plot the network
        net.3dplot(allnets[[i]][[j]][[k]], ROIs = comms)
        
        # Change the view
        view3d(userMatrix = view, zoom = 0.7)
        
        # Get the visit number
        vis.num = as.numeric(strsplit(names(allnets)[i], "Visit")[[1]][2])
        
        # Save rgl plot
        rgl.snapshot(paste0("allnetworks/RSnet_vis", vis.num, "_sub", j, "_net", k, ".png"))
        
        # Close the rgl window
        rgl.close()
      }
    }
  }
}