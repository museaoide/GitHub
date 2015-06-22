approxcs <- function(x,y,xout) {
  ## linear interpolation.  x and xout must be real, y can be complex
  ## returns interpolated y evaluated at the points in xout
  stopifnot(is.numeric(x), is.numeric(xout)) 
  xs <- sort(x,index.return=TRUE)
  x <- xs$x
  y <- y[xs$ix]
  xos <- sort(xout)
  N <- length(xout)
  M <- length(x)
  if (is.complex(y)) yout <- complex(length(xos)) else yout <- numeric(0)
  j <- 1
  k <- 1
  while (xos[j]<x[1] && j <= N ) { yout[j] <- y[1]; j <- j+1 }
  for (k in 2:M ) {
    while ( xos[j] < x[k] && j <= N ) {
      yout[j] <- (y[k] - y[k-1]) * (xos[j] - x[k-1]) / (x[k] - x[k-1]) + y[k-1]; j <- j+1
    }
  }
  while ( j<=N ) { yout[j] <- y[M]; j <- j+1}
  ## xos is returned because it might have been changed by sorting.
  return( list(y=yout,x=xos) )
}
