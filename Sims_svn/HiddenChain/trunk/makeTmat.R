makeTmat <- function(pij) {
  ## In this version, pij are just the first nState-1 rows of the transMat
  ## matrix, and this function just returns the full matrix.  More generally
  ## pij could be any set of parameters from which this function would
  ## construct transMat.  The mapping from pij to transMat will affect
  ## not only this routine, but also pijDraw()
  nState <- (1 + sqrt(1 + 4 * length(pij)))/2
  transMat <- matrix(pij,nState-1,nState)
  transMat <- rbind(transMat,1 - matrix(1, 1, nState - 1) %*% transMat)
  return(transMat)
}
