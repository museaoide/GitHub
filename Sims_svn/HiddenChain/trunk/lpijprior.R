lpijprior <- function(alpha, pij) {
  ## this version is for the case of a pure Dirichlet prior on the elements of pij
  ## Must provide both the log prior density and the initial state probabilities
  pij <- makeTMat(alpha)
  nstate <- dim(pij)[1]
  lpijNorm <- sum(lgamma(alpha)) - apply(lgamma(alpha), FUN=sum, MAR=2)
  eps <- .Machine$double.eps
  if (any( pij < 2*eps && alpha > 1 + 2*eps)) return( -1e20 )
  ##------------ assuming we treat phat0 as the ergodic distn -----------
  pbar <- qr.solve(rbind(pij - diag(nstate), rep(1,nstate)), c(rep(0,nstate),1))
  lpdf <- sum( (alpha - 1) * log(pij)) - lpijNorm
  ## ---------------(could return a constant pbar instead, e.g.)--------------------
  return(list(lpdf=lpdf, phat0=pbar))
}
