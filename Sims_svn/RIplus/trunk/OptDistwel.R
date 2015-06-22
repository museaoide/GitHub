optDistwel <- function(x,F,g) {
  ## dm is degree of mu(y) polynomial exponent.
  ## g is marginal pdf for y
  ## F is the exp(lambda*U(thet,y)) matrix
  ## x[1:dm] are the coefficients on the
  ## orthogonal polynomials of degree 0:dm
  ## ---------------------------------------
  dh <- length(x)
  nh <- dim(F)[1]
  nm <- dim(F)[2]
  ph <- poly(1:nh,degree=dh)
  h <- as.vector( exp(ph %*% x))
  h <- h/sum(h)
  m <- as.vector(g/(h %*% F))
  umat <- matrix(0,nh,nm)
  umat[F>0] <- log(F[F>0])
  ## browser()
  eu <- h %*% (F*umat) %*% m
  if(is.finite(eu) && is.real(eu)){
    return(-eu)
  }else{
    return(1e50)
  }
}
  
