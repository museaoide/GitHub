pijDraw <- function(S, alpha, nState, lags) {
  ## S is of length equal to the sample.  We use only time=lag (one before sample start) onward.
  ## Transpose of this function's alpha is the alpha for rdirichlet.  It defines the Dirichlet prior for pij.
  ## if (!require(MCMCpack)) stop("Program requires rdirichlet() from the MCMCpack package") # replaced by rdirichletCS.R
  Scount <- matrix(0, nState, nState)
  nT <- length(S)
  s12 <- cbind(S[lags:(nT - 1)], S[(lags + 1):nT])
  for (i1 in 1:nState) {
    for (i2 in 1:nState) {
      mtch <- s12 == matrix(c(i1,i2),nT-lags,2,byrow=TRUE)
      Scount[i1,i2] <- sum( ( mtch  %*% c(1,1) )==2)
    }
  }
  ##for (it in (lags + 1) : nT) {
  ##  Scount[S[it - 1], S[it]] = Scount[S[it - 1], S[it]] + 1
  ##}
  pij <- t(rdirichletCS(nState, Scount + t(alpha)))
  return(pij)
}
