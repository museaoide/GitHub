ldnorm <- function( x, mean, rooti) {
  ## rooti is inverse of cholesky factor of variance matrix
  n <- length(x)
  xt <- (x - mean) %*% rooti
  llh <- -.5 * n * log(2*pi) + sum(log(diag(rooti))) - .5 * sum(xt^2)
  return(llh)
}
