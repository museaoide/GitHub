panelLH <- function(ydata, lags, By, SigFacs, cs, Lambda, Mu, OmFac) {
  ## ydata:    periods by variables by countries (groups)
  ## By:       VAR coefficients, equations by variables by lags
  ## SigFacs:  nc inverses of square roots of covariance matrices for country residuals: SigFac %*% Sigma %*% t(SigFac) = diag(nv)
  ##           nv by nv by nc array.  These are inverses of cholesky factors, so their determinant is the product of diagonals.
  ## cs:       nv by nc matrix of country dummy coefficients
  ## Lambda:   lags by nv by nv matrix of coefficients connecting country const vectors to means of country initial conditions vectors
  ## OmFac:    inverse of square root of covariance matrix of initial conditions, as lags*nv by lags*nv array.
  require(abind)                        #Once this is packaged, the "require()"'s on every call should be dropped
  require(tensor)
  T <- dim(ydata)[1]
  nv <- dim(ydata)[2]
  nc <- dim(ydata)[3]
##   Xmat <- array(0, c(T-lags, nc, nv, lags))
##   for(ix in (lags+1)1:T)
##     Xmat[ix-lags, , ] <- t(ydata[ix - (1:lags), , drop=FALSE])
  Bmat <- aperm(By,c(3,2,1))
  Bmat <- Bmat[seq(lags,1,by=-1),,]
  Bmatd1 <- lags * nv
  Bmatd2 <- nv
  dim(Bmat) <- c(Bmatd1, Bmatd2)
  resid <- array(0, c(T - lags, nv, nc))
  for (it in 1:(T-lags))
    resid[it, , ] <- ydata[lags + it, , ] - t(tensor(ydata[seq(it+lags-1, it, by=-1), ,], By, alongA=c(1,2), alongB=c(3,2))) - cs
  resid <- tensor(resid, SigFacs, alongA=2, alongB=2)  # ends up T-lags by nc by nv
  resid0 <- ydata[1:lags, , ] - tensor(Lambda, cs, alongA=3, alongB=1) - array(Mu, c(lags, nv, nc))
  resid0 <-  tensor(resid0, array(OmFac, c(lags,nv,lags,nv)), alongA=c(1,2), alongB=c(1,2)) # ends up lags by nc by lags by nv
  llh <- sum( log( apply(SigFacs, 3, diag))) + sum(diag(OmFac)) - .5* (sum( resid^2) + sum(resid0^2)) - .5*log(2*pi)*T*nv*nc
  return(llh)
}
