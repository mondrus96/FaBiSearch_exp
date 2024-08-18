########################################################################################
### This script contains functions for doing analysis for the simulation experiments ###
########################################################################################

##### Load the Hausdorff function #####
source("Hausdorff.R")

##### Load the simulations functions #####
source("../../sharedfunctions/GenerateFunctions.R")
source("../../sharedfunctions/Simulations.R")

##### 1) Load all the data and organize final results for FaBiSearch change point detection #####
FBScpts = function(sims = 1:8){
  # sims = which simulation numbers to evaluate
  
  # Define the output variable which will be a dataframe of the results
  output = list()
  # Loop through the different simulations
  for(i in sims){
    # List the output files
    simouts = list.files(paste0("../../fabisearch/cpd/Sim", i), pattern = paste0("^sim", i, ".*.rda$"))
    
    # Loop through the different simulation output files and append
    results = list()
    for(j in simouts){
      # Load the results
      load(paste0("../../fabisearch/cpd/Sim", i, "/", j))
      
      # Add to the results
      results = rbind(results, all.results)
      }
    # Define the current output
    curr.out = list()
    
    # Find the average rank
    mean.rank = mean(results[!duplicated(results$i), "rank"])
    curr.out$K = mean.rank
    
    # Find the average compute time
    mean.cmput.T = mean(results[!duplicated(results$i), "mins"])
    curr.out$computeT = mean.cmput.T
    
    # Find the final change points
    final.cpts = results[results$stat_test, c("i", "T")]
    colnames(final.cpts)[1] = "rep"
    rownames(final.cpts) = NULL
    curr.out$changepoints = final.cpts
    
    # Append the current output to the final output
    output[[i]] = curr.out
  }
  # Return the final results
  return(output)
}

##### 2) Load all the data and organize final results for NCPD change point detection #####
NCPDcpts = function(sims = 1:8, nsims = 100){
  # sims = which simulation numbers to evaluate
  # nsims = number of simulation repetitions to be evaluated
  
  # Define the output variable which will be a dataframe of the results
  output = list()
  # Loop through the different simulations
  for(i in sims){
    # Read the output files, if they exist
    if (file.exists(paste0("../../NCPD/cpd/Sim", i, "/sim", i, "SignificantCPs.txt"))){
      final.cpts = read.table(paste0("../../NCPD/cpd/Sim", i, "/sim", i, "SignificantCPs.txt"))
      final.cpts = final.cpts[,1]
      
      # Calculate the difference between change points - if 
      # negative, then this is the start of a new simulation
      simnumdiff = diff(final.cpts) < 0
      # Create a vector of simulation numbers based on this pattern
      simnum = 1
      for(j in 2:length(simnumdiff)){
        simnum = c(simnum, ifelse(simnumdiff[j], simnum[length(simnum)] + 1, simnum[length(simnum)]))
      }
      # Need to add to beginning because of using diff()
      simnum = c(simnum[1], simnum)
      # Append this to the currcpts results
      final.cpts = cbind(simnum, final.cpts)
      
      # Rename columns
      colnames(final.cpts) = c("rep", "T")
    } else {
      # Initialize an empty matrix
      final.cpts = matrix(NA, 1, 2)
      
      # Rename columns
      colnames(final.cpts) = c("rep", "T")
      
      # Remove the row
      final.cpts = final.cpts[-1,]
    }
    
    # Load the extra data
    load(paste0("../../NCPD/cpd/Sim", i, "/sim", i, "ExtraResults.rda"))
    
    # Define the current output
    curr.out = list()
    
    # Find the average cluster size used
    curr.out$K = extraResults$tryK
    
    # Find the average compute time
    curr.out$computeT = extraResults$computeT/nsims
    
    # Add the final change points, and add as dataframe
    rownames(final.cpts) = NULL
    curr.out$changepoints = as.data.frame(final.cpts)
    
    # Append the current output to the final output
    output[[i]] = curr.out
  }
  # Return the output
  return(output)
}

##### 3) Calculate final results for all simulations #####
computeResults = function(resultslist, type, dists = c(1, 10), Haus = 10, nsims = 100){
  # resultslist = output of the FBScpts or NCPDcpts function, contains change point results
  # type = "FBS" for FaBiSearch, "NCPD" for NCPD
  # dists = distances to judge for TP/FP calculations, is a vector of numeric values
  # Haus = distance to use as maximum for Hausdorff calculations
  # nsims = number of sims
  
  # Initialize the output variable
  output = c()
  # Loop through the results list
  for(i in 1:length(resultslist)){
    # Loop through the different distance metrics
    acc.output = c()
    for(j in dists){
      # Define the boundaries of TP
      bounds = mapply(seq, (SIMtruecpts[[i]] - j), (SIMtruecpts[[i]] + j))
      if (i %in% c(7, 8)){
        bounds = as.matrix(unique(as.vector(bounds)))
      }
      
      # If change points exist in the results, then calculate the TP/FP rate(s)
      TPrate = FPrate = c()
      if(length(resultslist[[i]]$changepoints) > 0){
        # If there are change points, then continue to calculate TP rate
        if(!is.null(ncol(bounds))){
          # Loop through each column, or each change point
          for(k in 1:ncol(bounds)){
            # Calculate the TP rate, and keep value at maximum of 1 - using the min() function
            TPrate = c(TPrate, min(1, sum(resultslist[[i]]$changepoints$T %in% bounds[, k])/nsims))
          }
        } else {
          # Assign TPrate 0 if there are no change points
          TPrate = NA
        }
        FPrate = c(FPrate, sum(!resultslist[[i]]$changepoints$T %in% as.vector(bounds))/nsims)
      } else {
        # If no change points were detected by the algorithm(s),
        # then just assign a TP and FP rate of 0
        TPrate = FPrate = 0
      }
      
      # Append columns together and add to the output
      if(length(TPrate) > 1){
        FPrate = c(FPrate, rep("", length(TPrate) - 1))
      }
      currout = cbind(TPrate, FPrate)
      colnames(currout) = c(paste0("TP", j), paste0("FP", j))
      acc.output = cbind(acc.output, currout)
    }
    # If there are no change points, then trueT is NA
    if(is.null(SIMtruecpts[[i]])){
      trueT = NA
    } else {
      trueT = SIMtruecpts[[i]]
      if(all(diff(trueT) == 1) & length(trueT) > 1){
        trueT = paste0(min(trueT), ":", max(trueT))
      }
    }
    
    # Define bounds for Hausdorff distance and initialize the output Hdist
    bounds = mapply(seq, (SIMtruecpts[[i]] - Haus), (SIMtruecpts[[i]] + Haus))
    Hdist = c()
    # If there are change points, then continue to calculate Hausdorff distance
    if(length(bounds) != 0){
      # Loop through each simulation rep, or each change point
      for(k in unique(resultslist[[i]]$changepoints$rep)){
        Hcpts = resultslist[[i]]$changepoints$T[resultslist[[i]]$changepoints$rep == k]
        Hcpts = Hcpts[Hcpts %in% bounds]
        Hdist = c(Hdist, Hausdorff(Hcpts, SIMtruecpts[[i]], SIMlengths[i]))
      }
    } else {
      # Assign Hdist of NA if there are no change points
      Hdist = NA
    }
    # Take the mean
    Hdist = round(mean(Hdist), 4)
    
    # Append the other relevant information
    outlen = nrow(acc.output)
    # If the number of change points is greater than 1, need to append blank values
    # for all of the remaining results
    Sim = i
    MeanK = round(resultslist[[i]]$K, 2)
    computemins = round(as.numeric(resultslist[[i]]$computeT), 2)
    if(outlen > 1){
      Hdist = c(Hdist, rep("", outlen - 1))
      computemins = c(computemins, rep("", outlen - 1))
      MeanK = c(MeanK, rep("", outlen - 1))
    }
    
    sim.output = cbind(Sim = Sim, MeanK = MeanK, trueT = trueT, acc.output, 
                       Hdist = Hdist, computemins = computemins)
    
    # Append to final results
    output = rbind(output, sim.output)
  }
  # Return the output
  output = as.data.frame(output)
  
  # Call the saveTables function
  saveTables(output, type, dists)
  return(output)
}

##### 4) Save tables from the computeResults function #####
saveTables = function(output, type, dists = c(1, 10)){
  # output = output from the computeResults function
  # type = the type of change point detection method used - either "FBS" or "NCPD"
  # dists = distances to judge for TP/FP calculations, is a vector of numeric values
  
  # Change simulation numbering so that it is not repeated
  for(i in output$Sim){
    outlen = sum(output$Sim == i)
    if(outlen > 1){
      output$Sim[output$Sim == i] = c(i, rep("", outlen - 1))
    }
  }
  
  # Prepare outputs for LaTeX
  if (type == "FBS"){
    # Change column names for LaTeX
    accnames = c()
    for(a in 1:length(dists)){
      accnames = c(accnames, paste0("TP ", dists[a]), paste0("FP ", dists[a]))
    }
    colnames(output) = c("Sim.", "mean rank", "t*", accnames, "Haus. dist.", "Compute mins.")
    # Save it as a .csv file
    write.csv(output, "FBS.csv", row.names = FALSE)
  } else if (type == "NCPD"){
    # Change column names for LaTeX
    accnames = c()
    for(a in 1:length(dists)){
      accnames = c(accnames, paste0("TP ", dists[a]), paste0("FP ", dists[a]))
    }
    colnames(output) = c("Sim.", "K", "t*", accnames, "Haus. dist.", "Compute mins.")
    # Save it as a .csv file
    write.csv(output, "NCPD.csv", row.names = FALSE)
  }
}


##### 4) Evaluating the stationary network estimates for FaBiSearch #####
netResults = function(sims = 1:8, FBS, SIMlengths){
  # sims = which simulation numbers to evaluate
  # FBS = list of FaBiSearch results, from the FBScpts() results
  # SIMlengths = vector of simulation lengths
  
  # Define the output variable which will be a dataframe of the results
  output = data.frame("Accuracy")
  colnames(output) = "Sim."
  # Loop through the different simulations
  for(i in sims){
    # Print the current simulation being evaluated
    print(paste0("Simulation: ", i))
    
    # Define breaks the indices of stationary segments
    if(i == 1){
      SIMbounds = list(1:200)
    } else if (i == 2){
      SIMbounds = list(1:100, 101:200)
    } else if (i == 3){
      SIMbounds = list(1:100, 101:200, 201:300, 301:400)
    } else if (i == 4){
      SIMbounds = list(1:200, 201:400, 401:600)
    } else if (i == 5){
      SIMbounds = list(1:100, 101:200, 201:300)
    } else if (i == 6){
      SIMbounds = list(1:100, 101:200, 201:300)
    } else if (i == 7){
      SIMbounds = list(1:100, 151:200)
    } else if (i == 8){
      SIMbounds = list(1:100, 106:200)
    }
    
    # Retrieve the stationary networks and rename to simnets object
    load(paste0("../../fabisearch/statnets/Sim", i, "/sim", i, "_nets.rda"))
    simnets = all.results
    rm(all.results)
    
    # Define the simacc variable which will store the current simulations accuracy results
    simacc = c()
    
    # Loop through the individual simulations
    for(j in 1:length(simnets)){
      # Define the boundaries for stationary segments
      currbounds = c(0, FBS[[i]]$changepoints$T[FBS[[i]]$changepoints$rep == j],
                     SIMlengths[i])
      
      # Print the current simulation repetition being evaluated
      print(paste0(j))
      
      # Retrieve the current simulation sigmas, apply cutoffs
      SIMsigmas = do.call(paste0("sim", i), list(j, "sigma"))
      SIMsigmas = lapply(SIMsigmas, ">", 0.5)
      SIMsigmas = lapply(SIMsigmas, "*", 1)
      
      # Loop through the individual stationary networks
      for(k in 1:length(simnets[[j]])){
        # Define the current network
        currnet = simnets[[j]][[k]]
        currnet = currnet[upper.tri(currnet)]
        
        # Define the indices of the current segment
        currseg = (currbounds[k] + 1):currbounds[k + 1]
        
        # Find the greatest area of overlap based on time points
        overlap = c()
        for(l in 1:length(SIMbounds)){
          overlap[l] = sum(currseg %in% SIMbounds[[l]])
        }
        currSIMsigma = SIMsigmas[[match(max(overlap), overlap)]]
        currSIMsigma = currSIMsigma[upper.tri(currSIMsigma)]
        
        # Compare the upper triangle between currSIMsigma and the currnet
        simacc = c(simacc, sum(currnet == currSIMsigma)/length(currnet))
      }
    }
    # Append the results to the output
    output = cbind(output, round(mean(simacc)*100, 2))
    colnames(output)[i + 1] = i
    
    # Save the results as a .csv file for LaTeX
    write.csv(output, file = "graphingresults.csv", row.names = FALSE)
  }
}

##### 5) Creating the csv file for the comparison between FaBiSearch and NCPD #####
plotComp = function(FBS, NCPD, sims = 1:8, dist = 10){
  # FBS = dataframe of FaBiSearch results
  # NCPD = dataframe of NCPD results
  # sims = simulations to include in the output
  # dist = maximum distance used for output, coincides with TPrate/FPrate
  
  # Merge results into one dataframe
  allresults = rbind(cbind(Method = "FaBiSearch", FBS), cbind(Method = "NCPD", NCPD))
  
  # Select only the relevant columns
  allresults = allresults[, c("trueT", "Method", "Sim", paste0("TP", dist), paste0("FP", dist))]
  
  # Reorganize the dataframe
  output = c()
  for(i in sims){
    # Start with empty vectors
    FBS.res = NCPD.res = RateType = c()
    
    # Add the TP rate(s), if applicable
    if(!is.null(SIMtruecpts[[i]])){
      FBS.res = c(FBS.res, as.numeric(allresults$TP10[allresults$Method 
                                                      == "FaBiSearch" & allresults$Sim == i]))
      NCPD.res = c(NCPD.res, as.numeric(allresults$TP10[allresults$Method 
                                                        == "NCPD" & allresults$Sim == i]))
      RateType = c(RateType, allresults$trueT[allresults$Sim == i & 
                                                allresults$Method == allresults$Method[1]])
    }
    
    # Define the FP rate(s) and the grouping column
    FBS.res = c(FBS.res, as.numeric(allresults$FP10[allresults$Method 
                                                    == "FaBiSearch" & allresults$Sim == i][1]))
    NCPD.res = c(NCPD.res, as.numeric(allresults$FP10[allresults$Method 
                                                      == "NCPD" & allresults$Sim == i][1]))
    RateType = c(RateType, "FP")
    
    # Bind these separate columns together and add into output
    output = rbind(output, cbind(Sim = i, RateType, NCPD.res, FBS.res))
  }
  
  # Save the results as a csv file
  write.csv(output, "resultsforplot.csv", row.names = FALSE)
}