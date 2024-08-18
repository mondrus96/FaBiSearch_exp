# Load functions and data
source("tb_functions.R")
load("views.rda")
load("story.rda")

### Organize and save results ###
# Set working directory to the results folder
setwd("../results")
# Load the change point and network data
allnets = tbnets(1:8, 1:4)
allcpts = tbcpts(1:8, 1:4)
# Look at change point density
cptdens = cptdens(allcpts)
# Find maximums and print
cptmax = c()
for(i in 1:4){
  cptmax = c(cptmax, cptdens[[i]]$x[cptdens[[i]]$y
                                    == max(cptdens[[i]]$y)])
}
# Find the second maximum in the third run
r3dens = cptdens[[3]]
cptmax = c(cptmax, 
           r3dens$x[r3dens$x < 850][max(r3dens$y[r3dens$x < 850]) == 
                                      r3dens$y[r3dens$x < 850]])
cptmax = round(cptmax[order(cptmax)])
print(cptmax)
for(i in 1:length(cptmax)){
  print(i)
  print(story[(cptmax[i]*4 - 25):(cptmax[i]*4 + 25)])
}
# Save change points for LaTeX file
HPcpt = allcpts[,c(1, 4)]
rownames(HPcpt) = NULL
# Save as a csv file
write.csv(HPcpt, file = "HPcpt.csv", row.names = FALSE)
rm(HPcpt)

### Look at aggregate networks ###
# Look at overall/aggregate brain networks
aggnets(allnets, views)

### Look at degree distributions ###
# Estimate degree distributions
degdists = getdegdist(allcpts, allnets)
# Plot and save degree distributions for run 3
for(i in 1:8){
  plotdegdist(degdists, i, 3)
}

### Look at inter and intrasubject similarities ###
# Get similar inter and intra subject networks
similnets = findsimilnets(allnets)
# Keep only intrasubject comparisons
intrasubnets = similnets[similnets$sub1 == similnets$sub2,]
print(intrasubnets[1:20,])
# Keep only intersubject comparisons
intersubnets = similnets[similnets$sub1 != similnets$sub2,]
print(intersubnets[1:20,])
# Plot and save all networks
plotall(allnets, views[[1]])