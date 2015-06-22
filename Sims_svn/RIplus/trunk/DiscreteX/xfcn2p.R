xfcn2p <- function(ygivenx, y) {
  a <- .5  #cross-effect of other price on demand
  require(tensor)
  ybar <- ygivenx %*% y
  nx <- dim(ygivenx)[1]
  ny <- dim(y)[1]
  my <- dim(y)[2]
  x <- .5 *ybar + 1/(2 * (1 - a))
  return(x)
}
