panelPrior <- function (yinit, xinit=NULL, makexdum=TRUE, lags, lambda, mnprior=list(tight=.2,decay=.5),vprior=list(sig=1,w=1)) {
  ## yinit:     Tinit x nvar x ngroup array of initial conditions for y  (Usually Tinit=lags)
  ## xinit:     Tinit x nx x ngroup array of initial country x's.  (maybe just country dummies)
  ## lambda:    overall scale on the co-persistence dummies
  ## sig:       like vprior.sig in rfvar3:  approximate prior modes for std deviations of disturbances
  require(abind)
  nvar <- dim(yinit)[2]
  ngroup <- dim(yinit)[3]
  if (makexdum) {
    xinit <- abind(array(rep(c(diag(ngroup)), each =lags), c(lags, ngroup, ngroup)), xinit, along=3)
  }
  nx <- dim(xinit)[2]
  ybar <- apply(yinit, c(2,3), mean)
  xbar <- apply(xinit, c(2,3), mean)
  ydum <- array(0, c(lags+1, nvar, ngroup))
  xdum <- array(0, c(lags+1, nx, ngroup))
  for (il in 1:(lags+1)) {
    ydum[il, , ] <- ybar
    xdum[il, , ] <- xbar
  }
  for (iv in 1:nvar) {
    ydum[ , iv, ] <- lambda * ydum[ , iv, ] / vprior$sig[iv]
  }
  pbreaks <- (lags+1) * 1:(ngroup-1)
  ydum <- aperm(ydum, c(1,3,2)) # lags, group, var
  ydum <- matrix(ydum, ncol=nvar)
  xdum <- aperm(xdum, c(1,3,2))
  xdum <- matrix(xdum, ncol=nx)
  varp <- varprior(nv=nvar, nx=nx, lags=lags, mnprior=mnprior, vprior=vprior)
  pbreaks <- c(pbreaks, varp$pbreaks)
  ydum <- rbind(ydum, varp$ydum)
  xdum <- rbind(ydum, varp$ydum)
  return(list( ydum=ydum, xdum=xdum, pbreaks=pbreaks))
}
