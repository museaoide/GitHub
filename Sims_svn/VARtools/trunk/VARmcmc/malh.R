malh <- function(parvec,xdata) {
  ##  np <- length(parvec)-1
  ##  T <- length(xdata)
  ##  a <- parvec[1 : np]
  ##  mu <- parvec[np+1]
  ##  rx <- filter(c(rep(0,np-1),rev(a),rep(0,np-1)),a,sides=1)
  ##  rx <- rx[(2*np-1):(3*np-2)]
  ##  sig <- toeplitz(c(rx,rep(0,T-np)))
  ##  sch <- t(chol(sig))
  ##  six <- solve(sch,xdata)
  ##  onex <- solve(sch,rep(1,T))
  ##  muhat <- sum(six * onex)/sum(onex^2)
  ##  xtil <- xdata-muhat
  ##  ldetsch <- sum(log(diag(sch)))
  ##  sixtil <- six - muhat * onex
  ##  llh <- -.5 * T * log(2 * pi) -ldetsch -.5 * sum(sixtil^2) - .5 * (mu-muhat)^2 * sum(onex^2)
  ## ------------  version above treats mu as free; version below computes llh concentrated wrt mu
  np <- length(parvec)
  a <- parvec
  T <- length(xdata)
  ## rx <- filter(c(rep(0,np-1),rev(a),rep(0,np-1)),a,sides=1)
  rx <- convolve(a,a,type="open")
  rx <- rx[np : (2 * np - 1)]
  sig <- toeplitz(c(rx,rep(0,T-np)))
  sch <- t(chol(sig))
  six <- solve(sch,xdata)
  onex <- solve(sch,rep(1,T))
  muhat <- sum(six * onex)/sum(onex^2)
  xtil <- xdata-muhat
  ldetsch <- sum(log(diag(sch)))
  sixtil <- six - muhat * onex
  llh <- -.5 * T * log(2 * pi) -ldetsch -.5 * sum(sixtil^2)
  cat(sum(sixtil^2),xtil %*% solve(sig,xtil),"\n")
  return(c(-llh,muhat))
}
  
