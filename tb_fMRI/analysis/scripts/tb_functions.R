################################################################################################
##### This script contains functions for doing analysis for the task-based fMRI experiment #####
################################################################################################

# Load libraries
library(reshape2)
library(fabisearch)
library(rgl)
library(ggplot2)
library(ggridges)
library(extrafont)
library(dplyr)

# Change base font of ggplot plots
theme_set(theme_gray(base_size = 20, base_family = "CMU Serif"))

##### 1) Load all the data and organize final networks #####
tbnets = function(subjects, runs){
  # subjects = which subjects to estimate networks for
  # runs = which runs to estimate networks for
  
  # Define the final nets variable which will store the finalized networks
  final.nets = list()
  # Loop through the defined visits
  for(i in subjects){
    # Print the current subject
    print(paste0("Subject: ", i))
    
    # Load results, based on visit number
    load(paste0("../../statnets/Sub", i, "/HPnets", i, ".rda"))
    all.nets = all.results
    rm(all.results)
    
    # Loop through 
    subj.nets = list()
    for(j in runs){
      # Print the current subject
      print(paste0("Run: ", j))
      
      # Initialize the run nets list
      run.nets = list()
      
      # Sub-select the current run
      curr.run = all.nets[[j]]
      
      # Loop through the stationary blocks
      for(k in 1:length(curr.run)){
        curr.stat = curr.run[[k]]
        
        # Define the iterative loop variables for the while loop
        num.edges = 0
        l = length(curr.stat)
        while(num.edges < 100){
          # Select the current network
          curr.net = curr.stat[[l]]
          
          # Calculate and save the number of edges
          num.edges = sum(curr.net[lower.tri(curr.net)])
          
          # Reduce l by one
          l = l - 1
        }
        run.nets[[k]] = curr.net
      }
      subj.nets[[j]] = run.nets
    }
    final.nets[[i]] = subj.nets
  }
  # Return the final nets object
  return(final.nets)
}

##### 2) Load all the change points into a data frame #####
tbcpts = function(subjects, runs){
  # subjects = which subjects to estimate networks for
  # runs = which runs to estimate networks for

  # Loop through the defined visits
  HPall = c()
  for(i in subjects){
    # Load results, name based on visit number
    load(paste0("../../cpd/Sub", i, "/HPcpd", i, ".rda"))
    
    # Add a column to denote the subject
    all.results$subject = i
    
    # Append to the HPall object, and remove the all.results object
    HPall = rbind(HPall, all.results)
    rm(all.results)
  }
  # Keep only final change points
  HPall = HPall[HPall$stat_test,]
  
  # Keep only necessary columns and rename
  HPall = HPall[,c(5,1:3)]
  colnames(HPall)[c(2,4)] = c("run", "time") 
  
  # Return the results
  return(HPall)
}

##### 3) Look at overall/aggregate brain network across time points #####
aggnets = function(final.nets, views){
  # final.nets = list of finalized nets
  
  # Compare the node connections across subjects/in general
  sum.adj.matr = matrix(0, 333, 333)
  for(i in 1:length(final.nets)){
    for(j in 1:length(final.nets[[i]])){
      for(k in 1:length(final.nets[[i]][[j]])){
        sum.adj.matr = sum.adj.matr + final.nets[[i]][[j]][[k]] 
      }
    }
  }
  
  ### Heatmap plotting ###
  # Open png file
  png("agg_heatmap.png", width = 1000, height = 1000, res = 250)
  # Display a heatmap
  heatmap(sum.adj.matr, Rowv = NA, Colv = NA, cexRow = 0.5, cexCol = 0.5)
  # Close png file
  dev.off()
  
  # Print the nodes with the greatest number of edges, using average node degree
  node.deg = colSums(sum.adj.matr)/length(final.nets)
  print("5 Highest degree nodes (descending order)")
  print(names(node.deg[order(-node.deg)][1:5]))
  
  # Print the percent of nodes localized to the right side
  print("Percentage of nodes laterized to the right,")
  cutoffs = c(5, 10, 20)
  for(i in cutoffs){
    # Printe the degree cutoff
    print(paste0("degree >", i))
    
    # Define the current proportion
    curr.nodes = names(node.deg[node.deg > i])
    curr.nodes = as.numeric(sapply(strsplit(curr.nodes, "V"), "[[", 2))
    print(sum(curr.nodes >= 162)/length(curr.nodes))
  }
  
  ### 3D network plotting ###
  # Melt the matrix, and order by the number of connections
  melt.sum = sum.adj.matr
  melt.sum[upper.tri(melt.sum)] = NA
  melt.sum = melt(melt.sum, na.rm = TRUE)
  print(melt.sum[order(-melt.sum$value)[1:10],])
  
  # Look at the quantiles and print
  edge.quant = quantile(melt.sum$value, probs = 0.999)
  print(edge.quant)
  
  # Plot based on quantiles
  plot.matr = 1*(sum.adj.matr >= edge.quant)
  net.3dplot(plot.matr)
  
  # Loop through the different views
  for(i in 1:length(views)){
    # Change the view
    view3d(userMatrix = views[[i]], zoom = 0.7)
    
    # Save rgl plot
    rgl.snapshot(paste0("agg_net", i, ".png"))
  }
  
  # Close the rgl window
  rgl.close()
}

##### 4) Plot a histogram of density of change points #####
cptdens = function(allcpts){
  # allcpts = all change point data
  
  # Load fonts
  loadfonts(device = "win")
  
  # Define the binsize as well as cutoffs for grouping of runs
  bin.size = 25
  cutoffs = list(seq(-50, 350, bin.size),
                 seq(300, 750, bin.size),
                 seq(680, 980, bin.size),
                 seq(900, 1400, bin.size))
  
  # Initialize the change point density estimate outputs
  cptdensout = list()
  
  # Open the png file
  png(filename = "hp_cps_density.png", width = 2000, height = 1300, res = 200)
  # Create a window of plots
  par(mfrow=c(1,4), oma = c(4.5, 5.5, 3, 3), mar=c(1,1,1,1))
  for(i in 1:4){
    run = i
    dens = density(allcpts$time[allcpts$run == run])
    datalen = length(allcpts$time[allcpts$run == run])
    h = hist(allcpts$time[allcpts$run == run], breaks = cutoffs[[run]], prob = FALSE, ylim = c(0, 8),
         xlab = "", ylab = "", main = "", cex.axis = 1.5, family = "CMU Serif")
    lines(dens$x, datalen*(1/sum(h$density))*dens$y, type = "l", col = "red")
    title(paste0("Run ", i), cex.main = 2, family = "CMU Serif")
    
    # Append the densities to the cptdensout list
    cptdensout[[i]] = dens
  }
  title(xlab = "Time", ylab = "Frequency", outer = TRUE, cex.lab = 3, family = "CMU Serif")
  # Close the png file
  dev.off()
  
  # Return the change point density estimates
  return(cptdensout)
}

##### 5) Plot degree distributions over time #####
getdegdist = function(allcpts, allnets){
  # allcpts = all change points
  # allnets = all stationary networks
  
  # Loop through
  all.results = c()
  for(i in 1:length(allnets)){
    # Loop through all subjects
    for(j in 1:length(allnets[[i]])){
      # Add temporal boundaries to change points depending on run
      if(j == 1){
        curr.cpts = c(0, allcpts[allcpts$subject == i & allcpts$run == j, "time"], 324)
      } else if (j == 2){
        curr.cpts = c(325, allcpts[allcpts$subject == i & allcpts$run == j, "time"], 661)
      } else if (j == 3){
        curr.cpts = c(662, allcpts[allcpts$subject == i & allcpts$run == j, "time"], 925)
      } else if (j == 4){
        curr.cpts = c(926, allcpts[allcpts$subject == i & allcpts$run == j, "time"], 1290)
      }
      curr.nets = allnets[[i]][[j]]
      
      # For the current collection of networks, loop through
      for(k in 1:length(curr.nets)){
        # Define the midpoint
        y = ceiling(mean(curr.cpts[k:(k+1)]))
        x = 0:20
        
        # Find the degree of all nodes
        v.degs = colSums(as.matrix(curr.nets[[k]]))
        height = c()
        for(l in 0:20){
          prop = sum(v.degs == l)/333
          height[l+1] = 1-exp(-10*prop)
        }
        
        # Append the midpoint, subject, and run data
        curr.results = cbind(i, j, x, y, height)
        
        # Append to the overall matrix
        all.results = rbind(all.results, curr.results)
      }
    }
  }
  all.results = as.data.frame(all.results)
  rownames(all.results) = NULL
  colnames(all.results)[1:2] = c("sub", "run")
  
  # Return final result
  return(all.results)
}

##### 6) Plot the degree distribution #####
plotdegdist = function(degdists, sub, run){
  # degdists = output from the degdist() function
  # sub = subject number
  # run = run number
  
  # Define temporal boundaries based on run
  if (run == 1){
    r.min = 0
    r.max = 350
  } else if (run == 2){
    r.min = 350
    r.max = 700
  } else if (run == 3){
    r.min = 700
    r.max = 950
  } else if (run == 4){
    r.min = 950
    r.max = 1300
  }
  
  # Select the distributions to plot
  curr.selection = degdists[degdists$sub == sub & degdists$run == run, 3:5]
  
  # Plot the selected degree distributions
  ggplot(curr.selection, aes(x, y, height = height*30, group = y)) + ylim(r.min, r.max) + xlab("Node Degree") + ylab("Time") +
    geom_ridgeline(fill = "lightblue") + theme(panel.background = element_blank(), text = element_text(size = 20),
                                               axis.line = element_line(colour = "black"))
  # Save as a png
  ggsave(filename = paste0("ridgeplots/degdist_sub", sub, "_run", run, ".png"), width = 10, height = 8)
}

##### 7) Finding similarities with and between subjects #####
findsimilnets = function(allnets){
  # allnets = all stationary networks
  
  # Initialize the output
  matr.simil = c()
  # Loop through the subjects
  for(i in 1:length(allnets)){
    # Select the current subjects networks
    sub.i = allnets[[i]]
    
    # Loop through all subjects
    for(j in 1:length(allnets)){
      sub.j = allnets[[j]]
      for(k in 1:4){
        sub.i.k = sub.i[[k]]
        sub.j.k = sub.j[[k]]
        for(l in 1:length(sub.i.k)){
          for(m in 1:length(sub.j.k)){
            curr.net.i = sub.i.k[[l]]
            curr.net.j = sub.j.k[[m]]
            
            # Calculate the similarity
            net.simil = sum(curr.net.i*curr.net.j)/2
            
            # Append results to dataframe, and order based on which subject is lower in value
            if (i < j){
              matr.simil = rbind(matr.simil, cbind(k, i, l, j, m, net.simil))
            } else if (j < i){
              matr.simil = rbind(matr.simil, cbind(k, j, m, i, l, net.simil))
            } else if (i == j){
              # If the two subjects are the same - as in this is intra subject, order based on block
              if (m < l){
                matr.simil = rbind(matr.simil, cbind(k, j, m, i, l, net.simil))
              } else if (l < m){
                matr.simil = rbind(matr.simil, cbind(k, i, l, j, m, net.simil))
              }
            }
          }
        }
      }
    }
  }
  # Rename the columns
  colnames(matr.simil)[1:5] = c("run", "sub1", "block1", "sub2", "block2")
  
  # Remove duplicate values
  matr.simil = distinct(as.data.frame(matr.simil))
  
  # Order descending by network similarity and remove rownames
  matr.simil = matr.simil[order(-matr.simil[, "net.simil"]),]
  rownames(matr.simil) = NULL
  
  # Return the similarity matrix
  return(matr.simil)
}

##### 8) Plot and save all networks #####
plotall = function(allnets, view){
  # allnets = all stationary networks
  # view = view from which to plot the rgl window from
  
  # Loop through subjects
  for(i in 1:length(allnets)){
    # Loop through runs
    for(j in 1:length(allnets[[i]])){
      # Loop through networks
      for(k in 1:length(allnets[[i]][[j]])){
        # Plot the network
        net.3dplot(allnets[[i]][[j]][[k]])
    
        # Change the view
        view3d(userMatrix = view, zoom = 0.7)
        
        # Save rgl plot
        rgl.snapshot(paste0("allnetworks/HPnet_sub", i, "_run", j, "_net", k, ".png"))
        
        # Close the rgl window
        rgl.close()
      }
    }
  }
}