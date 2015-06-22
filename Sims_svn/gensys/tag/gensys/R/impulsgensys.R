#' impulsgensys
#'
#' Impulse responses for a solved linear RE model
#'
#' @param gout a list of system matrices as produced by \code{gensys}
#' @param horiz number of periods over which to calculate responses
#' @param shocksize if non-NULL, a list of sizes of initial shocks
#'
#' @return a 3-d array, with dimensions variables x shocks x timespan, with
#'    the first two dimensions named.
#' @seealso \code{gensys}
#' @export
#' 
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
