HCdraw <- function(psGYfilt=matrix(0, 5, 2), transMat=diag(2), lags) {
  # draws states conditional on full data set.
  ## psGYfilt:    filtered state probabilities, a nT x nState matrix, p(S_t | Y_t)
  ## transMat:    Markov transition matrix.  Columns add to one.
  ## lags:        number of y lags
  nT <- dim(psGYfilt)[1]
  nState <- dim(psGYfilt)[2]
  states <- rep(0,nT)
  p <- psGYfilt[nT, ]
  for (it in seq(nT, lags, -1)) {
    pd <- cumsum(c(0, p))
    states[it] <- findInterval(runif(1), pd, all.inside=TRUE)
    p <- psGYfilt[it-1, ] * transMat[states[it], ]
    psum <-sum(p)
    if (psum > 0) {
      p <- p/sum(p)
    } else {
      p <- rep(1/nState, nState)
      warning(paste("All states zero prob in HCdraw at", it))
    }
  }
  lPdf <- sum(log(transMat[cbind(states[(lags + 1):nT], states[lags:(nT - 1)])])) # likelihood of the drawn states
  return(list(states = states, lPdf = lPdf))
}
