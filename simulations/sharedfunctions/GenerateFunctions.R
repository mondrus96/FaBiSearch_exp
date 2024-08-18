# Auxiliary functions in simulation

# libraries required
library("MASS")

exp.decay = function(rho, p){
	temp = matrix(0, p, p)
	for(i in 1:p){
		for(j in i:p){
			temp[i,j] = temp[j, i] = rho^{abs(i-j)}
		}
	}
	return(temp)
}

new.mod = function(i, j){
	temp = i %% j
	if(temp > 0) return(temp)
	if(temp == 0) return(j)
}

group.corr = function(small.n, rho, p){
	temp = matrix(0, p, p)
	for(i in 1:p){
		for(j in i:min((i+ (small.n - new.mod(i, small.n))), p)){
			temp[i, j] = temp[j, i] = (i==j)*1 + (i!=j)*rho
		}
	}
	return(temp)
}

