#' Cosine transform
#' 
#' Creates a cosine transform matrix of the given size.
#' 
#' The result F is a real orthonormal matrix \code{crossprod(F)=crossprod(t(F) = diag(n))} with first row a constant
#' and succesive rows \code{k} after that oscillating  with periods \code{2n/(k-1)}.  .
#'
#' @param n Dimension of the matrix.
ctmat <- function(n) {
   F <- matrix(0,n,n)
  for (j in 1:n) for (k in 1:n) F[j, k] <- cos(  pi * (k -.5) * (j-1) / n)
  F <- c(1/sqrt(2), rep(1, n-1)) * F
  return(F)
}