blkDglz <- function(A, div=c(.99, 1-100*sqrt(.Machine$double.eps)), ctOrder=FALSE) {
  ## Returns P, D such that D is block diagonal (and triangular within
  ## blocks), the abs(roots) of the upper block of D are all less than div[1],
  ## the abs(roots) of the lower block are all greater than div(2), and the
  ## abs(roots) of the middle block are in between.  P %*% D %*% solve(P) = A.
  if ( !is.loaded("zhseqr")) dyn.load("/usr/lib/liblapack.so")
  sf <- schur(A)
  ## sf <- zgeesWrap(A)
  Q <- sf$Q
  T <- sf$T
  n <- dim(Q)[1]
  ev <- diag(as.matrix(sf$T))
  if (ctOrder) {
    comp <- function(a){Re(a) < div[1]}
  } else {
    comp <- function(a){abs(a) < div[1]}
  }
  nlow <- sum(comp(ev))
  sfwork <- schdiv(sf, comp=comp)
  ## ulndx <- comp(diag(as.matrix(sf$T)))
  ulndx <- if(nlow > 0) 1:nlow else NULL
  P <- matrix(sfwork$Q[ , ulndx], n, nlow)
  if (ctOrder) {
    comp <- function(a) { Re(a) >= div[1] & Re(a) <= div[2] }
  } else {
    comp <- function(a) { abs(a) >= div[1] & abs(a) <= div[2] }
  }
  nmid <- sum(comp(ev))
  ulndx <- if( nmid > 0 ) 1:nmid else NULL
  sfwork <- schdiv(sf, comp=comp)
  ## ulndx <- comp(diag(as.matrix(sf$T)))
  P <- cbind(P, sfwork$Q[ , ulndx])
  if (ctOrder) {
    comp <- function(a) { Re(a) > div[2] }
  } else {
    comp <- function(a) { abs(a) > div[2] }
  }
  nhi <- sum( comp(ev) )
  ulndx <- if (nhi > 0) 1:nhi else NULL
  sfwork <- schdiv(sf, comp=comp)
  ## ulndx <- comp(diag(as.matrix(sf$T)))
  P <- cbind(P, sfwork$Q[ , ulndx])
  Pinv <- solve(P)
  D <- Pinv %*% A %*% P
  return(list(P=P, Pinv=Pinv, D=D, blockDims=c(nlow, nmid, nhi)))
}
