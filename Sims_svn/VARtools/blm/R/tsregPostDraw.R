tsregPostDraw <- function(blmout, n) {
  cf <- blmout$lsout$coefficients
  ncf <- length(cf)
  s2draw <- with(blmout, postdf/rgamma(n, postdf))
  bdraw <- matrix(rnorm(n * ncf), ncol=ncf)
  bdraw <- s2draw * bdraw %*% chol(blmout$vcv) + matrix(cf, n, ncf, byrow=TRUE)
  s2draw <- with(blmout, s2draw * postssr/(2 * postdf))
  return(list(bdraw = bdraw, s2draw = s2draw))
}