rsf2csf <- function(sf) {
  T <- sf$T
  Q <- sf$Q
  n <- dim(T)[1]
  if (n==1) return(sf)
  ld <- diag(T[2:n, 1:(n-1), drop=FALSE])
  cd <- diag(T)
  icx <- which(abs(ld) > sqrt(abs(cd[1:(n-1)] * cd[2:n])) * 100*.Machine$double.eps)
  ## icx is the indexes of elements of the first diagonal below the main
  ## diagonal that are non-zero
  for (ild in icx+1) {
    W <- rsf2csf22(T[(ild - 1):ild, (ild - 1):ild])
    T[(ild - 1):ild, ] <- W %*% T[(ild-1):ild, ]
    T[ , (ild-1):ild] <- T[ , (ild -1):ild] %*% t(Conj(W))
    Q[ , (ild-1):ild] <- Q[ , (ild-1):ild] %*% t(Conj(W))
  }
  return(list(Q=Q, T=T))
}
rsf2csf22 <- function(A) {
  a <- A[1,1]
  b <- A[1,2]
  c <- A[2,1]
  d <- A[2,2]
  w2bar <- (a - d + sqrt(as.complex((a - d)^2 + 4 * b * c)))/(2 * b)
  sf <- sqrt(1 + abs(w2bar)^2)
  w1 <- 1/sf
  w2bar <- w2bar/sf
  W <- matrix(c(w1, -w2bar, Conj(w2bar), w1),2)
  return(W)
}
