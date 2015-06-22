HClhElement <- function(y, pyGs=function(y=0, S=c(1, 2), ...){}, phat, transMat=diag(2), ...) {
  ## pyGS:  log pdf of y(t) given S(t) and past data and parameters given in ...  Must return a vector of the same dimension as its S argument
  ##        (which can be scalar).
  ## transMat:  Markov transition matrix.  p(t+1)=transMat %*% p(t)
  ## ...:   exactly those variables needed as extra arguments to pyGs
  nState <- dim(transMat)[2]
  plh <- pyGs(y, S=1:nState, ...)
  phatnew <- exp(plh) * (transMat %*% phat)
  element <- log(sum(phatnew))
  phatnew <- phatnew / sum(phatnew)
  if (any(is.nan(phatnew))) phatnew <- rep(1,nState)/nState
  return(list(lhElement=element,phatnew=phatnew))
}
