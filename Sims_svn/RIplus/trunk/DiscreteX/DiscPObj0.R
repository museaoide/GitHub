DiscPObj <- function(param, gy, y, U, alph) {
  ## param:  First nx-1 are marginal probabilities on x points,
  ##         remainder are the positions of the discrete x points.
  ## gy:     The marginal probabilities on y
  ## y:      The discrete y points
  ## U:      A function of x and y defining utility
  ## alph:   One over the utility cost of information
  ## obj:    The objective function value.
  ## "p" and "x" index the variables wrt which we are taking derivs
  ## "weight" indexes the nx rows of the joint density
  ## "y" indexes the ny columns of the joint density
  require("tensor")
  nx <- (length(param)+1)/2
  ny <- length(gy)
  p <- param[1:(nx - 1)]
  p <- matrix(c(p, 1 - sum(p)), 1,nx, dimnames=list(NULL, p=NULL))
  x <- matrix(param[nx:(2 * nx - 1)], 1, nx, dimnames=list(NULL, x=NULL))
  if ( any(p < 0) ) return(1e20)
  mm <- function(z) { matrix(z, nx, ny, byrow=TRUE) }
  sy <- function(z) {apply( z, MAR=1, FUN=sum)}
  Umat <- outer(c(x), c(y), FUN=U)
  ## peaU <- p %*% exp(alph * Umat)
  eaU <- exp(alph * Umat)
  dimnames(eaU) <- list(weight=NULL, y=NULL)
  DPpeaU <- eaU            
  names(dimnames(DPpeaU))[1] <- "p"
  peaU <- tensor(p, eaU, 2,1)
  DXU <- matrix(attr(Umat,"gradient"), nx, ny, dimnames=list(x=NULL, y=NULL))
  ##assumes U fcn written to return gradient attribute,
  ##e.g. via deriv().
  DXpeaU <- mm(tensor(p, alph * eaU, 2, 1)) * DXU
  h <-  gy / peaU
  DPh <- -DPpeaU * mm(c(gy / peaU^2))
  DXh <- -DXpeaU * mm(c(gy / peaU^2))
  eaUh <- eaU * mm(c(h))
  f <- c(p) * eaUh
  DPf <- array(0, c(nx,nx,ny), dimnames=list(p=NULL, weight=NULL, y=NULL))
  for (ip in 1:nx){
    for (iw in 1:nx) {
      DPf[ip, iw, ] <- p[iw]*eaU[iw, ] * DPh[ip, ] 
    }
    DPf[ip,ip, ] <- DPf[ip,ip, ] + eaUh[ip, ]
  }
  DXf <- array(0,c(nx,nx,ny), dimnames=list(x=NULL, weight=NULL, y=NULL))
  for (ix in 1:nx) {
    for (iw in 1:nx) {
      DXf[ix,iw, ] <- p[iw] * eaU[iw, ] * DXh[ix, ]
    }
    DXf[ix, ix, ] <- DXf[ix, ix, ] + p[ix] * eaUh[ix, ] * alph * DXU[ix, ] * h
  }
  DPobj <- tensor(DPf, Umat, c(2,3), c(1,2))
  DXobj <- tensor(DXf, Umat, c(2,3), c(1,2)) + sy(f * DXU * alph)
  ## Above are derivatives of the expected utility part.  On to the info cost part.
  pnew <- sy(f)
  DPpnew <- tensor(DPf, rep(1, ny), 3, 1)
  DXpnew <- tensor(DXf, rep(1, ny), 3, 1)
  ipplus <- p > 0
  fplus <- f[ipplus, ]
  pplus <- pnew[ipplus]
  roweight <- sy(eaUh)
  ygivenx <- c(1/roweight) * eaUh
  DProweight <- tensor(DPh, eaU, 2, 2)  #p=nx x w=nx matrix
  DXroweight <- diag(sy(eaU * alph * DXU)) + tensor(DXh, eaU, 2, 2)
  DPygivenx <- array(0,c(nx,nx,ny))
  DXygivenx <- array(0,c(nx,nx,ny))
  ## result of line below is weight=nx by p=nx by y=ny tensor.  Weight indexes rows of ygivenx, p indexes what deriv's are wrt.
  for (iw in 1:nx) {
    DPygivenx[ , iw, ] <- DProweight[ , iw] %o% (-(1/roweight[iw]^2) * eaUh[iw, ])
  }
  for (iy in 1:ny) {
    DPygivenx[ , , iy] <- DPh[ , iy] %o% (c(1/roweight) * eaU[ , iy])
  }
  for (iw in 1:nx) {
    DXygivenx[ , iw, ] <- DXroweight[ , iw] %o% ((-1/roweight[iw]^2) * eaUh[iw, ])
    DXygivenx[iw, iw, ] <- DXygivenx[iw,iw,] + (1/roweight[iw]^2) * eaUh[iw, ] * DXU[iw, ] *alph 
  }
  for (iy in 1:ny) {
    DXygivenx[ , , iy] <- DXygivenx[ , , iy] + DXh[ , iy] %o% (c(1/roweight) * eaU[ , iy])
  }                                                           
  obj <- sum(f  * Umat ) - (1/alph) * pnew %*% sy(log(ygivenx) * ygivenx)
  DPobj <- c(DPobj) - c((1/alph) * DPpnew %*% sy(log(ygivenx) * ygivenx)) + c(tensor(DPygivenx, c(pnew) * (1 + log(ygivenx)), 2:3, 1:2))
  DXobj <- c(DXobj) - c((1/alph) * DXpnew %*% sy(log(ygivenx) * ygivenx)) + c(tensor(DPygivenx, c(pnew) * (1 + log(ygivenx)), 2:3, 1:2))
  obj <- -obj                           #as input to minimizer
  attr(obj, "gradient") <- -c(cbind(diag(nx - 1), rep(-1, nx - 1)) %*% DPobj, DXobj) #recongizing that there are only nx-1 p parameters, with
                                        # n'th being 1 minus the rest.
  ## browser()
  ##attr(obj, "pnew") <- pnew             #allows checking that pnew = p (as it should at optimum)
  return(obj)
}
