impulsgensys <- function(gout, horiz, shocksize=NULL) {
  ## gout is output list from gensys
  if (!is.null(shocksize)) impact <- gout$impact * diag(shocksize) else impact <- gout$impact
  nv <- dim(gout$impact)[1]
  nshock <- dim(gout$impact)[2]
  impulse <- array(0, c(nv, nshock, horiz))
  impulse[ , , 1] <- impact
  for ( it in 2:horiz) {
    impulse[ , , it] <- gout$G1 %*% impulse[ , , it-1]
  }
  dimnames(impulse) <- c(dimnames(gout$impact), list(t=c(0:(horiz-1))))
  return(impulse)
}
