pyGsAR1 <- function(y, S, yl, x=NULL, parvec) { # rho, gam, sigsq, x=NULL) {
  ## y:      y(t)
  ## S:      states at which to evaluate (vector 1:nState, usually)
  ## yl:     lagged y's to use
  ## rho:    AR coefficients
  ## sigsq:  vector of sigsq values for the different states
  ## x:      constant and/or exogenous variables
  sigsq <- parvec$sigsq
  llh <- -.5 * (y - parvec$rho %*% yl - (if (!is.null(x)) x %*% parvec$gam else 0))^2 / sigsq[S] - .5 * log(sigsq[S])
  return(llh)
}
