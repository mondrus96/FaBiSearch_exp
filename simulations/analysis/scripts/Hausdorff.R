Hausdorff <- function(a, b, n_s){
  l_a = l_b = 0
  if(length(a)==0){
    a = c(0)
  }else{
    if(length(b)==0){
      b = c(0)
    }
  }
  for (i in 1:length(a)) {
    l_a_i = abs(a[i] - b)
    l_a = max(l_a, min(l_a_i))
  }
  for (i in 1:length(b)) {
    l_b_i = abs(b[i] - a)
    l_b = max(l_b, min(l_b_i))
  }
  return(max(l_a,l_b)/n_s)
}

