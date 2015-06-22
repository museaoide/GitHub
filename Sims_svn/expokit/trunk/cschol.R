cschol <- function(x) {
  ## like chol, but 60 times slower.  The only advantage is that it
  ## produces u.t. w with crossprod(w) = x even when x is p.s.d. but not p.d.
  ## Code by Christopher Sims, this version 7/22/2003.
  ##------------------------------------------------
  eps <- .Machine$double.eps
  n0 <- dim(x)[1]
  m <- dim(x)[2]
  if (!n0 == m) error("non-square argument")
  nerr <- 0
  for (i in 1:n0) {
    if (i > 1) x[i,1:(i-1)] <- 0
    w <- x[i:n0,i:n0, drop=FALSE]
    er <- 0
    n <- n0-i+1
    if (w[1,1] > 1000 * eps) {
      w[1,1] <- Re(sqrt(w[1,1]))
      if (n >= 2) {
        w[1,2:n] <- w[1,2:n]/w[1,1]
        wr <- w[1,2:n, drop=FALSE]
        w[2:n,2:n] <- w[2:n,2:n] - crossprod(wr)
      }
    } else {
      er <- 1
      if (w[1,1]<0) {
        if (w[1,1]<-1000*eps) {
          er <- -1
        }
        w[1,1] <- 0
      }
      if (n >= 2) w[2:n,1] <- 0
    }
    x[i:n0,i:n0] <- w
    if (er) {
      if (nerr==0) nerr <- i*er
    }
  }
  w <- x
  attr(w,"nerr") <- nerr
  return(w)
}
