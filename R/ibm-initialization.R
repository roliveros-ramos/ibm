# Initialization ----------------------------------------------------------

pop.initial = function(N, n, min=0, max=1) {
  out = matrix(runif(n*N,min=min,max=max), ncol=n, nrow=N)
  return(out)
}
