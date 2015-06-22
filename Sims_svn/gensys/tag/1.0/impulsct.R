impulsct <- function(impact, G1, interval, span) {
  ## impact and G1 are assumed to be of the form in the output of gensysct.R
  ## interval is the time interval at which impulse responses are to be evaluated.
  ## span is the time span over which the impulse responses are to be evaluated.
  ## to get the forward part of a gensysct solution, use -fmat as G1, fwt as impact,
  ## and then premultiply by ywt (i.e. tensor(ywt, resp, 2, 1) or array( ywt %*% matrix(resp, nrow=nv), c(nv,ns,nresp))
  nv <- dim(G1)[1]
  ns <- dim(impact)[2]
  nrsp <- ceiling(span/interval)+1
  resp <- array(0,c(nv, ns, nrsp))
  A <- padm(G1 * interval)
  resp[ , , 1] <- impact
  for (it in 2:nrsp) resp[ , , it] <- A %*% resp[ , , it-1]
  dimnames(resp) <- list(dimnames(impact)[[1]], dimnames(impact)[[2]], NULL)
  return(resp)
}
  
