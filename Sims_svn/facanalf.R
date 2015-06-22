facanalf <- function(d, n, s) {
  ## version to be argument to csminwel (returning only llh)
  ## d: vector of log std deviations of idiosyncratic components
  ## n: number of factors
  ## s: covariance matrix to be analyzed
  ## llh: 1/T times log of likelhood
  ## F: matrix of factor loadings: E[s]=diag(d) %*% (F %*% F' + I) %*% diag(d)
  d <- exp(d)
  sd <- s / (d %o% d)
  nv <- length(d)
  umv <- svd(sd)
  if (umv$d[n] < 1){
    llh <- -1e20
    F <- matrix(0,nv,n)
  } else {
    llh <- -sum(log(d)) - .5 * sum(log(umv$d[1:n])) - .5 * (sum(diag(s) / d^2) - sum(umv$d[1:n]) + n)
    F <- umv$v[,1:n] %*% diag(sqrt(umv$d[1:n]-1))
  }
  if (!is.null(dimnames(s))) dimnames(F) <- list(dimnames(s)[[1]],NULL)
  ## return(list(llh=llh,F=F,m=umv$d))  # for outside csminwel
  return(-llh)
}
