# Project y onto L2 ball of radius c
proc = function(y, c){
  y = as.numeric(y)
  if (sqrt(sum(y^2)) <= c) return(y)
  c * y / sqrt(sum(y^2))
}