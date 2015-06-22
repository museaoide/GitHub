arlh <- function(parvec,xdata) {
  T <- length(xdata)
  rho <- 2*atan(parvec[1])/pi
  sigsq <- exp(parvec[2])
  ## rx <- filter(c(rep(0,np-1),rev(a),rep(0,np-1)),a,sides=1)
  rx <- rho^(0:(T-1))*sigsq
  sig <- toeplitz(rx)
  sch <- t(chol(sig))
  six <- solve(sch,xdata)
  onex <- solve(sch,rep(1,T))
  muhat <- sum(six * onex)/sum(onex^2)
  ## xtil <- xdata-muhat
  ldetsch <- sum(log(diag(sch)))
  sixtil <- six - muhat * onex
  llh <- -.5 * T * log(2 * pi) -ldetsch -.5 * sum(sixtil^2)
  return(-llh)
}
