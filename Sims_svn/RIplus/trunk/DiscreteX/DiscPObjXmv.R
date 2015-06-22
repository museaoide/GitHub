DiscPObjXmv <- function(param, gy, y, U, alph, ...) {
  ## DiscPObj, but with y and x both potentially multivariate
  ## 
  ## param:  First is nx.  Next nx-1 are marginal probabilities on 1st n-1 x points,
  ##         remainder fill the nx x mx matrix of positions of the discrete x points,
  ## gy:     The ny marginal probabilities on y
  ## y:      The discrete y points, as an ny by my matrix
  ## U:      A function of x and y defining utility
  ## alph:   One over the utility cost of information
  ## obj:    The objective function value.
  ##
  require("tensor")
  ny <- dim(y)[1]
  my <- dim(y)[2]
  nx <- param[1]
  mx <- (length(param) - nx)/nx
  p <- param[2:nx]
  p <- matrix(c(p, 1 - sum(p)), 1, nx, dim=list(NULL, p=NULL))
  x <- matrix(param[nx + (1:(nx*mx))], nx, mx, dimnames=list(NULL, x=NULL))
  if ( any(p < 0) ) {
    print(paste("p < 0.  min", min(p)))
    p <- pmax(p,0)
  }
  mm <- function(z) { matrix(z, nx, ny, byrow=TRUE) }
  sy <- function(z) {apply( z, MAR=1, FUN=sum)}
  entel <- function(z) ifelse(z > 0, -z * log(z), 0)
  ## Umat <- outer(c(x), c(y), FUN=U, ...)
  Umat <- matrix(0, nx, ny)
  for (ix in 1:nx) {
    for (iy in 1:ny) {
      Umat[ix, iy] <- U(x[ix,], y[iy,], ...)
    }
  }
  eaU <- exp(alph * Umat)
  dimnames(eaU) <- list(weight=NULL, y=NULL)
  peaU <- tensor(p, eaU, 2,1)
  #if (any(peaU < 1e-290 * gy)) return(1e20)
  h <-  gy / peaU
  eaUh <- eaU * mm(c(h))
  f <- c(p) * eaUh
  pnew <- sy(f)
  ipplus <- pnew > 0
  ixp <- (1:nx)[ipplus]
  nr <- sum(ipplus)
  fplus <- f[ipplus, ,drop=FALSE]
  pplus <- pnew[ipplus]
  roweight <- sy(eaUh)[ipplus]
  ygivenx <- matrix(0, nx, ny)
  ygivenx[ipplus, ] <- c(1/roweight) * eaUh[ipplus, ,drop=FALSE]
  yent <- sum(entel(gy))
  obj <- sum(f  * Umat ) + (1/alph) * pnew %*% sy(entel(ygivenx))  - (1/alph) * yent 
  ## obj <- -obj                           #as input to minimizer
  info = -sum(entel(f)) + sum(entel(pnew)) + yent
  return(list(obj=obj, pnew=pnew, ygivenx=ygivenx, h=h, f=f, info=info, EU = sum(f * Umat)))
}
