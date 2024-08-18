# Each function takes in the iteration number (iter) which affects the seed number 
# and subsequently random number generation for the multivariate normal processes

# For each simulation, there are two inputs:

# iter = the iteration number - usually corresponds to the loop/repetition number
# return = by default is "y", the entire multivariate simulation. However, if set
#          to "sigma", will return a list of the original covariance matrix values

##### Simulation 1 #####
sim1 = function(iter, return = "y"){
  # Time points per segment
  t       = 200
  # Nodes per cluster
  s.p     = 200
  # Clusters
  K       = 2
  n.tot   = K * s.p
  # Within cluster correlation
  rho.wi  = 0.75
  # Between cluster correlation
  rho.bw  = 0.2
  # Covariance of first segment
  Sigma1 	      = group.corr(small.n = s.p, rho = rho.wi, p = n.tot)
  ind 	      = which(Sigma1 == 0)
  Sigma1[ind]   = rho.bw
  
  # Set up the simulation
  set.seed(iter*3)
  reorder = sample(c(1:n.tot), n.tot)
  Sigma2  = Sigma1[reorder, reorder]
  Y1      = mvrnorm(n = t, mu = rep(0, n.tot), Sigma = Sigma1)
  y 		  = Y1
  
  # Return the output
  if(return == "y"){
    # Return the simulation
    return(y)
  } else if(return == "sigma"){
    # Return original covariance matrices as list
    return(list(Sigma1))
  }
}
##### Simulation 1 #####

##### Simulation 2 #####
sim2 = function(iter, return = "y"){
  # Time points per segment
  t       = 100
  # Nodes per cluster
  s.p     = 200
  # Clusters
  K       = 2
  n.tot   = K * s.p
  # Within cluster correlation
  rho.wi  = 0.75
  # Between cluster correlation
  rho.bw  = 0.2
  # Covariance of first segment
  Sigma1 	      = group.corr(small.n = s.p, rho = rho.wi, p = n.tot)
  ind 	      = which(Sigma1 == 0)
  Sigma1[ind]   = rho.bw
  
  # Set up the simulation
  set.seed(iter*7)
  reorder = sample(c(1:n.tot), n.tot)
  Sigma2  = Sigma1[reorder, reorder]
  set.seed(iter*3)
  Y1      = mvrnorm(n = t, mu = rep(0, n.tot), Sigma = Sigma1)
  set.seed(iter*4)
  Y2 		  = mvrnorm(n = t, mu = rep(0, n.tot), Sigma = Sigma2)
  y 		  = rbind(Y1, Y2)
  
  # Return the output
  if(return == "y"){
    # Return the simulation
    return(y)
  } else if(return == "sigma"){
    # Return original covariance matrices as list
    return(list(Sigma1,
                Sigma2))
  }
}
##### Simulation 2 #####

##### Simulation 3 #####
sim3 = function(iter, return = "y"){
  # Inputs for time points 1-100
  s.p     = 200
  K       = 3
  n.tot   = K * s.p
  rho.wi  = 0.75
  rho.bw  = 0.2
  Sigma1  = group.corr(small.n = s.p, rho = rho.wi, p = n.tot)
  for(i in 1:n.tot){
    for(j in 1:n.tot){
      if(Sigma1[i,j] == 0){
        Sigma1[i,j] = rho.bw^(abs(i - j))
      }
    }
  }
  
  # Time points 101-200
  s.p     = 300
  K       = 2
  n.tot   = K * s.p
  rho.wi  = 0.75
  rho.bw  = 0.2
  Sigma2  = group.corr(small.n = s.p, rho = rho.wi, p = n.tot)
  for(i in 1:n.tot){
    for(j in 1:n.tot){
      if(Sigma2[i,j] == 0){
        Sigma2[i,j] = rho.bw^(abs(i - j))
      }
    }
  }
  
  # Set up the simulation
  t = 100
  set.seed(iter*2)
  Y1 = mvrnorm(n = t, mu = rep(0, n.tot), Sigma = Sigma1)
  set.seed(iter*11+100001)
  Y2 = mvrnorm(n = t, mu = rep(0, n.tot), Sigma = Sigma2)
  set.seed(iter*13)
  reorder = sample(c(1:n.tot), n.tot)
  Sigma3  = Sigma2[reorder, reorder]
  Y3 = mvrnorm(n = t, mu = rep(0, n.tot), Sigma = Sigma3)
  set.seed(iter*17+100001)
  reorder1 = sample(c(1:n.tot), n.tot)
  Sigma4  = Sigma1[reorder1, reorder1]
  Y4 = mvrnorm(n = t, mu = rep(0, n.tot), Sigma = Sigma4)
  y  = rbind(Y1, Y2, Y3, Y4)
  
  # Return the output
  if(return == "y"){
    # Return the simulation
    return(y)
  } else if(return == "sigma"){
    # Return original covariance matrices as list
    return(list(Sigma1,
                Sigma2,
                Sigma3,
                Sigma4))
  }
}
##### Simulation 3 #####

##### Simulation 4 #####
sim4 = function(iter, return = "y"){
  # Inputs for time points 1-200
  s.p		    = 400
  K 		    = 2
  n.tot 	  = K * s.p
  rho.wi 	  = 0.75
  rho.bw 	  = 0.2
  Sigma1 	  = group.corr(small.n = s.p, rho = rho.wi, p = n.tot)
  for(i in 1:n.tot){
    for(j in 1:n.tot){
      if(Sigma1[i,j] == 0){
        Sigma1[i,j] = rho.bw^(abs(i - j))
      }
    }
  }
  
  # 201-400
  reorder1 	= c((1:200), (401:600), (201:400), (601:800))
  Sigma2    = Sigma1[reorder1, reorder1]
  
  # 401-600
  reorder2 	= c((1:100), (401:500), (101:200), (501:600), (201:300), (601:700), (301:400), (701:800))
  Sigma3    = Sigma1[reorder2, reorder2]
  
  # Set up the simulation
  t = 200
  set.seed(iter*11)
  Y1 = mvrnorm(n = t, mu = rep(0, n.tot), Sigma = Sigma1)
  set.seed(iter*11+100001)
  Y2 = mvrnorm(n = t, mu = rep(0, n.tot), Sigma = Sigma2)
  set.seed(iter*13)
  Y3 = mvrnorm(n = t, mu = rep(0, n.tot), Sigma = Sigma3)
  y  = rbind(Y1, Y2, Y3)
  
  # Return the output
  if(return == "y"){
    # Return the simulation
    return(y)
  } else if(return == "sigma"){
    # Return original covariance matrices as list
    return(list(Sigma1,
                Sigma2,
                Sigma3))
  }
}
##### Simulation 4 #####

##### Simulation 5 #####
sim5 = function(iter, return = "y"){
  # Inputs for time segments
  t         = 100
  s.p		    = 100
  K 		    = 2
  n.tot 	  = K * s.p
  rho.wi 	  = 0.75
  rho.bw 	  = 0.2
  Sigma1 	  = group.corr(small.n = s.p, rho = rho.wi, p = n.tot)
  for(i in 1:n.tot){
    for(j in 1:n.tot){
      if(Sigma1[i,j] == 0){
        Sigma1[i,j] = rho.bw^(abs(i - j))
      }
    }
  }
  
  # Set up the simulation
  set.seed(iter*12345)
  reorder = sample(c(1:n.tot), n.tot)
  Sigma2  = Sigma1[reorder, reorder]
  set.seed(iter*11)
  Y1      = mvrnorm(n = t, mu = rep(0, n.tot), Sigma = Sigma1)
  set.seed(iter*4)
  Y2 		  = mvrnorm(n = t, mu = rep(0, n.tot), Sigma = Sigma2)
  set.seed(iter*3)
  Y3      = mvrnorm(n = t, mu = rep(0, n.tot), Sigma = Sigma1)
  y 		  = rbind(Y1, Y2, Y3)
  y  = rbind(Y1, Y2, Y3)
  
  # Return the output
  if(return == "y"){
    # Return the simulation
    return(y)
  } else if(return == "sigma"){
    # Return original covariance matrices as list
    return(list(Sigma1,
                Sigma2,
                Sigma1))
  }
}
##### Simulation 5 #####

##### Simulation 6 #####
sim6 = function(iter, return = "y"){
  # Inputs for time points 1-100, and 201-300
  t         = 100
  s.p		    = 100
  K 		    = 7
  n.tot 	  = K * s.p
  rho.wi 	  = 0.75
  rho.bw 	  = 0.2
  Sigma1 	  = group.corr(small.n = s.p, rho = rho.wi, p = n.tot)
  for(i in 1:n.tot){
    for(j in 1:n.tot){
      if(Sigma1[i,j] == 0){
        Sigma1[i,j] = rho.bw^(abs(i - j))
      }
    }
  }
  
  # Inputs for time points 101-200
  s.p		    = 350
  K 		    = 2
  n.tot 	  = K * s.p
  rho.wi 	  = 0.8
  rho.bw 	  = 0
  Sigma2 	  = group.corr(small.n = s.p, rho = rho.wi, p = n.tot)
  for(i in 1:n.tot){
    for(j in 1:n.tot){
      if(Sigma2[i,j] == 0){
        Sigma2[i,j] = rho.bw^(abs(i - j))
      }
    }
  }
  
  # Set up the simulation
  set.seed(iter*11)
  Y1      = mvrnorm(n = t, mu = rep(0, n.tot), Sigma = Sigma1)
  set.seed(iter*4)
  Y2 		  = mvrnorm(n = t, mu = rep(0, n.tot), Sigma = Sigma2)
  set.seed(iter*3)
  Y3      = mvrnorm(n = t, mu = rep(0, n.tot), Sigma = Sigma1)
  y 		  = rbind(Y1, Y2, Y3)
  y  = rbind(Y1, Y2, Y3)
  
  # Return the output
  if(return == "y"){
    # Return the simulation
    return(y)
  } else if(return == "sigma"){
    # Return original covariance matrices as list
    return(list(Sigma1,
                Sigma2,
                Sigma1))
  }
}
##### Simulation 6 #####

##### Simulation 7 #####
sim7 = function(iter, return = "y"){
  # Time points per segment
  t       = 150
  # Nodes per cluster
  s.p     = 200
  # Clusters
  K       = 2
  n.tot   = K * s.p
  # Within cluster correlation
  rho.wi  = 0.75
  # Between cluster correlation
  rho.bw  = 0.2
  # Covariance of first segment
  Sigma1 	      = group.corr(small.n = s.p, rho = rho.wi, p = n.tot)
  ind 	      = which(Sigma1 == 0)
  Sigma1[ind]   = rho.bw
  
  # Set up the simulation
  set.seed(iter*7)
  reorder = sample(c(1:n.tot), n.tot)
  Sigma2  = Sigma1[reorder, reorder]
  
  set.seed(iter*3)
  Y1      = mvrnorm(n = t, mu = rep(0, n.tot), Sigma = Sigma1)
  set.seed(iter*4)
  Y2 		  = mvrnorm(n = t, mu = rep(0, n.tot), Sigma = Sigma2)
  
  # Combine Y1 and Y2 over 50 time points smoothly
  Ymid = c()
  for(j in 1:50){
    weight = 1/(1+2.71828^(-0.2*(j-25)))
    Ymid = rbind(Ymid, Y1[100+j,]*(1-weight) + Y2[j,]*(weight))
  }
  
  # Combine into final
  y = rbind(Y1[1:100,], Ymid, Y2[51:150,])
  
  # Return the output
  if(return == "y"){
    # Return the simulation
    return(y)
  } else if(return == "sigma"){
    # Return original covariance matrices as list
    return(list(Sigma1,
                Sigma2))
  }
}
##### Simulation 7 #####

##### Simulation 8 #####
sim8 = function(iter, return = "y"){
  # Time points per segment
  t       = 150
  # Nodes per cluster
  s.p     = 200
  # Clusters
  K       = 2
  n.tot   = K * s.p
  # Within cluster correlation
  rho.wi  = 0.75
  # Between cluster correlation
  rho.bw  = 0.2
  # Covariance of first segment
  Sigma1 	      = group.corr(small.n = s.p, rho = rho.wi, p = n.tot)
  ind 	      = which(Sigma1 == 0)
  Sigma1[ind]   = rho.bw
  
  # Set up the simulation
  set.seed(iter*7)
  reorder = sample(c(1:n.tot), n.tot)
  Sigma2  = Sigma1[reorder, reorder]
  set.seed(iter*3)
  Y1      = mvrnorm(n = t, mu = rep(0, n.tot), Sigma = Sigma1)
  set.seed(iter*4)
  Y2 		  = mvrnorm(n = t, mu = rep(0, n.tot), Sigma = Sigma2)
  
  # Combine Y1 and Y2 over 50 time points smoothly
  Ymid = c()
  for(j in 1:5){
    Ymid = rbind(Ymid, Y1[100+j,]*(1-(j/5)) + Y2[j,]*(j/5))
  }
  
  # Combine into final
  y = rbind(Y1[1:100,], Ymid, Y2[51:150,])
  
  # Return the output
  if(return == "y"){
    # Return the simulation
    return(y)
  } else if(return == "sigma"){
    # Return original covariance matrices as list
    return(list(Sigma1,
                Sigma2))
  }
}
##### Simulation 8 #####