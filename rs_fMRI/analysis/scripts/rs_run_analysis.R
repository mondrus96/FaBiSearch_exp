# Load functions and data
source("rs_functions.R")
load("views.rda")

### Organize and save results ###
# Set working directory to the results folder
setwd("../results")
# Load the change point and network data
allnets = rsnets(2:3)
allcpts = rscpts(2:3)

### Look at inter and intrasubject similarities ###
# Plot heatmap for visit 2 of subject 1
vis2sub1 = allnets$Visit2[[1]]
for(i in 1:length(vis2sub1)){
  # Open a png
  png(paste0("vis2sub1_net", i, ".png"), width = 1000, height = 1000, res = 250)
  # Plot a heatmap
  heatmap(vis2sub1[[i]], col = grey(c(0.97,0.1)), 
          Rowv = NA, Colv = NA, labRow = FALSE, labCol = FALSE)
  # Close the png
  dev.off()
}
# Plot and save all networks
comms = c("Default")
plotall(allnets, views[[1]], comms)