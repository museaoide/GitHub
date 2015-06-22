HClh <- function(y, x=NULL, pyGs, phat0, transMat, lags, ...) {
  ## ... contains only arguments for pyGS that don't change with t.
  if ( is.null(dim(y)) ) {
    dim(y) <- c(length(y),1)
  }
  if (!is.null(x)) {
    if (is.null(dim(x)))
      dim(x) <- c(length(x),1)
    if (dim(x)[1] != dim(y)[1]) stop("x and y must have the same nrow")
  }
  nT <- dim(y)[1]
  psGYfilt <- matrix(NA, nT, length(phat0))
  psGYfilt[lags, ] <- phat0
  llhVec <- vector("numeric",nT)
  ## llh <- 0.0
  for (it in (lags+1):nT) {
    lhe <- HClhElement(y=y[it, ], pyGs=pyGs, phat=psGYfilt[it-1, ], transMat=transMat, yl=y[(it-1):(it-lags), ], x=x[it, ], ...)
    psGYfilt[it, ] <- lhe$phatnew
    ## llh <- llh + lhe$lhElement
    llhVec[it] <- lhe$lhElement
  }
  return(list(llhVec=llhVec, psGYfilt=psGYfilt))
}
