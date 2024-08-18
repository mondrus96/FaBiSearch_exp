# This script will run NCPD for simulations

# Get the start time
start.T = Sys.time()

# Loop through the 100 iterations of the simulation
for(i in 1:100){
  # Print the iteration
  print(i)
  
  # Set the seed
  set.seed(i*123)
  
  # Generate simulation depending on the sim type
  y = do.call(paste0("sim", simnum), list(i))
  
  # Inputs for change point function
  # Length of the input data/time course
  T = nrow(y)
  # Column matrix of time
  x                 = matrix(1:T,ncol=1)
  # Algorithm type (1 means the algorithm is exhausted first and then inference
  # is performed, 2 means the inference is performed after each cp found
  type              = 1          
  # First time point for the algorithm
  lower             = 0
  # Last time point for the algorithm
  upper             = T
  # List for output
  cp.index       = list()
  # Minimum distance between change points
  maxn              = 50
  # Quantile threshold to locate outliers
  quan              = 0.95
  # Number of bootstrap\permutation replicates
  n.boot           = 10
  # Size of average block length for stationary bootstrap
  block            = 0.2
  # Quantile for confidence bounds adjusting for multiple comparisons
  quantile1        = 0.01
  # Permutation/stationary bootstrap procedure (1 means permutation, 2 stationary
  # bootstrap
  permtype         = 2
  # Names for output (these can be changed)
  actualcps        = paste0("sim", simnum, "ActualCPs.txt")
  bootcps          = paste0("sim", simnum, "BootstrapCPs.txt")
  signcps          = paste0("sim", simnum, "SignificantCPs.txt")
  elapsed.time     = paste0("sim", simnum, "ElapsedTime.txt")
  # Set up the cluster
  np = detectCores(logical = FALSE)
  cl = makeCluster(getOption("cl.cores",np))
  
  # Run NCPD
  spectralfunction(x,y,type,try.K,T,lower,upper,cp.index,maxn,quan,n.boot,block,quantile1,permtype,actualcps,bootcps,signcps,elapsed.time)
  
  # Stop the cluster
  stopCluster(cl)
}

# Get the end time
end.T = Sys.time()

# Get time difference and save this
extraResults = list(computeT = difftime(end.T, start.T, units = "mins"),
                    tryK = try.K)
save(extraResults, file = paste0("sim", simnum, "ExtraResults.rda"))
