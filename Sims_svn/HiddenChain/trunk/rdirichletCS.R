rdirichletCS <- function (n, alpha) {
  ## This actually allows matrix alpha, like the writeup of rdirichlet says, but not like the
  ## actual rdirichlet in MCMCpack.  The rows of alpha correspond to separate draws.
  if (is.null(dim(alpha))) {
    l <- length(alpha)
  } else {
    l <- dim(alpha)[2]
    stopifnot(n == dim(alpha)[1])
  }
  x <- matrix(rgamma(l * n, t(alpha)), nrow = l)
  sm <-  rep(1, l) %*% x
  return(t(x)/as.vector(sm))
}
