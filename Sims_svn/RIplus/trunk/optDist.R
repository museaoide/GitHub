optDist <- function(x,F,g,dh,dm) {
  ## dm is degree of mu(y) polynomial exponent, dh for h(thet).
  ## g is marginal pdf for y
  ## F is the exp(lambda*U(thet,y)) matrix
  ## x[1:dh] are the coefficients on the orthogonal
  ## polynomials of degree 1:dh
  ## x[(dh+1):(dh+dm+1)] are the coefficients on the
  ## orthogonal polynomials of degree 0:dm
  ## ---------------------------------------
  if(is.null(dim(x))) {
    nx <- 1
    dim(x) <- c(length(x),1)
  }else {
    nx <- dim(x)[2]
  }
  ch <- x[1:dh,,drop=FALSE]
  cm <- x[(dh+1):(dh+dm+1),,drop=FALSE]
  nh <- dim(F)[1]
  nm <- dim(F)[2]
  pm <- poly(1:nm,degree=dm)
  pm <- cbind(rep(1/sqrt(nm),nm),pm)
  m <- exp(pm %*% cm)
  ph <- poly(1:nh,degree=dh)
  h <- exp(ph %*% ch)
  h <-t(t(h)/apply(h,2,sum,drop=FALSE))
  eh <- 1-(F %*% m) #nh rows
  ##browser()
  ehs <- apply(eh,2,sum)
  dim(ehs) <- c(1,nx)
  eh <- cbind(t(ehs),crossprod(eh,ph))
  em <- g - m * crossprod(F,h) #nm rows
  em <- crossprod(em , pm)
  ##browser()
  if(all(is.finite(c(em,eh)))) {
    return(t(cbind(eh,em[,1:dm,drop=FALSE])))
  }else {
    return(matrix(1e20,dm+dh+1,nx))
  }
}
  
