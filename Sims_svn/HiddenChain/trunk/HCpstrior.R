HCpstrior <- function(parvec, y, x, pyGs, lags, parPrior, ...) {
  ## To keep this function generic, we do not deduce lags from parvec, but make it a separate
  ## argument.  Only the pij component of parvec is parsed out here.
  lpijp <- lpijprior(alpha=parvec$alpha, pij=parvec$pij, phat0)     #log joint prior on pij,phat0
  lpijpfac <- lpijp$lpdf
  phat0 <- lpijp$pbar
  lpYpfac <- lparprior(parvec$pY)       #log prior on pYgS parameters
  nT <- dim(ydata)[1]
  if ( is.null(dim(y)) ) {
    dim(y) <- c(length(y),1)
    if (!is.null(x) && dim(x)[1] != dim(y)[1]) stop("x and y must have the same nrow")
  }
  nT <- dim(y)[1]
  psGYfilt <- phat0
  llh <- 0.0
  for (it in (lags+1):nT) {
    lhe <- HClhElement(y=y[it, ], pyGs=pyGs, phat=psGYfilt[it-1, ], transMat=transMat, yl=y[(it-lags):(it-1), ], x=x[it, ], parvec=parvec, )
    psGYfilt[it, ] <- lhe$phatnew
    llh <- llh + lhe$lhElement
  }
  dObsOut <- setDumObs(parPrior, sigsq, y0, x0)
  ## need separate prior for each state.  In principle, both rho-gam and sigsq could change.
  ## Is this trying to be general VAR switching?  Seems so.
  return(llh-lpij)
}
