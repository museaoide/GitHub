DynInfoObj <- function(x, rho, nu, lambda, beta, w, K=100, trigger=.1, lean=TRUE) {
  ## calculates w' Sigma w - logdet Sigma + logdet(rho * Sigma * rho' + nu),
  ## while penalizing approaches to zero of any eigenvalues of 
  n <- length(w)
  ## if (!is.null(dim(x))) {
  ##   nc <- dim(x)[2]
  ## } else {
  ##   nc <- 1
  ##   x <- matrix(x, n, nc)
  ## }
  deltaf <- matrix(0,n,n)
  ## objmat <- matrix(0, n, nc)
  ## for (ic in 1:nc) {
  ## deltaf[upper.tri(deltaf, diag=TRUE)] <- x[ , ic]
  deltaf[upper.tri(deltaf, diag=TRUE)] <- x
  delta <- crossprod(deltaf)
  sigma <- doubling(rho, nu - delta)
  cnd <- chol(nu - delta, pivot=TRUE, LINPACK=FALSE)
  ##this relies on chol with pivot returning negative numbers on diagonal for non-psd matrices.
  ## if (any(diag(cnd) <= 0)) {
  ## obj <- 1e20
  ## penalty <- trigger
  ## print("cnd trigger")
  ## browser()
  ## } else {
  omega <- rho %*% sigma %*% t(rho) + nu
  csigma <- chol(sigma, pivot=TRUE, FALSE)
  comega <- chol(omega, pivot=TRUE, FALSE)
  penalty1 <- diag(csigma)
  penalty1 <- penalty1[penalty1 < trigger]
  penalty2 <-abs( diag(deltaf))
  penalty2 <- penalty2[penalty2 < trigger]
  penalty <- c(penalty1, penalty2)
  if (is.null(penalty) ) {
    obj <- t(w) %*% sigma %*% w - lambda * sum(log(diag(csigma)))
    + lambda * beta * sum(log(diag(comega)))
    penalty <- trigger
  } else {
    if (any(penalty1 <= 0)) {
      obj <- 1e20
      penalty <- trigger
    } else {
      if(any(penalty2 <= 0)) {
        obj <- 1e20
        penalty <- trigger
      } else {
        obj <- t(w) %*% sigma %*% w - lambda * sum(log(diag(csigma)))
        + lambda * beta * sum(log(diag(comega)))
      }  
    }
    ##}
    if (lean ) return( obj + K * sum(((trigger - penalty)/penalty)^2))
    return(list(obj=obj, penalty=K * sum(((trigger - penalty)/penalty)^2), sigma=sigma,
                omega=omega, delta=delta))
  }
}
