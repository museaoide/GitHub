DiscCAPMobj <- function(param, gy, y, alph) {
  ## DiscPObj, but with useful extras returned, for use in analyzing results,
  ## not in csminwel.
  ## param:  First nx-1 are marginal probabilities on x points,
  ##         remainder are the positions of the discrete x points as (nq-1) x nx  array
  ##         The minus 1 reflects the need for x to sum to wealth down columns.
  ## gy:     The marginal probabilities on y, an ny x 1 array
  ## y:      The discrete y points, an ny x nq array
  ## --------- With nq assets, ny is likely to by prod(nyi, i=1:nq), so the y's are in an nq-d 
  ## --------- rectangular grid.  But this program doesn't need to know that.
  ## U:      A function of x and y defining utility
  ## alph:   One over the utility cost of information
  ## obj:    The objective function value.
  ## "p" and "x" index the variables wrt which we are taking derivs
  ## "weight" indexes the nx rows of the joint density
  ## "y" indexes the ny columns of the joint density
  ##------------------------
  ##
  require("tensorA")
  W <- 2                                #set wealth to 2, so it's easy to see that its components are not probabilities.
  nq <- dim(y)[2]
  ny <- dim(y)[1]
  SIGZ <- diag(c(rep(.04, nq-1), 0))             # covariance matrix of irreducible part.
  MUZ <- c(rep(.5, nq-1), .4)                    # mean of irreducible part
  nx <- (length(param) + 1) / (nq + 1)
  ny <- length(gy)
  p <- param[1:(nx - 1)]
  p <- matrix(c(p, 1 - sum(p)), 1,nx, dimnames=list(NULL, p=NULL))
  x <- matrix(param[nx:(nq * nx - 1)], nq-1, nx, dimnames=list(NULL, x=NULL))
  x <- rbind(x, W - apply(x, 2, sum))
  if ( any(p < 0) ) return(1e20)
  mm <- function(z) { matrix(z, nx, ny, byrow=TRUE) }
  sy <- function(z) {apply( z, MAR=1, FUN=sum)}
  entel <- function(z) ifelse(z > 0, z * log(z), 0)
  xy <- t(y %*% x)
  zmux <- matrix(MUZ %*% x, nx, ny)
  zsigxx <- apply(x * (SIGZ %*% x), 2, sum)
  zsigxx <- matrix(zsigxx, nx, ny)
  Umat <- xy + zmux - (xy + zmux)^2/2 -zsigxx/2
  yt <- rep(to.tensor(t(y), c(q=nq, y=ny)), nx, pos=3, names="x")
  zmuxt <- rep(to.tensor(zmux, c(x=nx, y=ny)), nq, pos=1, names="q")
  zsigxxt <- rep(to.tensor(zsigxx, c(x=nx, y=ny)), pos=1, names="q")
  DXU <- (yt + zmuxt) * (1 - rep(to.tensor(xy + zmux, c(x=nx, y=ny)), nq, names="q")) - rep(to.tensor(SIGZ %*% x, c(q=nq, x=nx)), ny, names="y")
  ## now DXU is an nq x nx x ny tensor
  ## Umat <- outer(c(x), c(y), FUN=U)
  ## peaU <- p %*% exp(alph * Umat)
  eaU <- exp(alph * Umat)
  dimnames(eaU) <- list(weight=NULL, y=NULL)
  DPpeaU <- eaU            
  names(dimnames(DPpeaU))[1] <- "p"
  ## peaU <- mul.tensor(p, eaU, 2,1)           # isn't this just p %*% eau ?
  peaU <-  p %*%  eaU
  if (any(peaU < 1e-290 * gy)) return(1e20)
  ## DXpeaU <- mm(tensor(p, alph * eaU, 2, 1)) * DXU
  DXpeaU <- c(p) * (alph * mul.tensor(to.tensor(eaU, c(x=nx, y=ny)), NULL, DXU, by = c("nx","ny")))
  ## working here 5/2/2012
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
    ## DXf[ix, ix, ] <- DXf[ix, ix, ] + p[ix] * eaUh[ix, ] * alph * DXU[ix, ] * h
    DXf[ix, ix, ] <- DXf[ix, ix, ] + p[ix] * eaUh[ix, ] * alph * DXU[ix, ] #the h on the end redundant (?)
  }
  DPobj <- tensor(DPf, Umat, c(2,3), c(1,2))
  DXobj <- tensor(DXf, Umat, c(2,3), c(1,2)) + sy(f * DXU )
  ## Above are derivatives of the expected utility part.  On to the info cost part.
  pnew <- sy(f)
  DPpnew <- tensor(DPf, rep(1, ny), 3, 1)
  DXpnew <- tensor(DXf, rep(1, ny), 3, 1)
  ipplus <- pnew > 0
  ixp <- (1:nx)[ipplus]
  nr <- sum(ipplus)
  fplus <- f[ipplus, ,drop=FALSE]
  pplus <- pnew[ipplus]
  roweight <- sy(eaUh)[ipplus]
  ygivenx <- c(1/roweight) * eaUh[ipplus, ,drop=FALSE]
  DProweight <- tensor(DPh, eaU[ipplus, ,drop=FALSE], 2, 2) #p=nx x w=nx matrix
  DXroweight <- diag(sy(eaUh * alph * DXU))[, ipplus,drop=FALSE] + tensor(DXh, eaU[ipplus, ,drop=FALSE], 2, 2)
  DPygivenx <- array(0,c(nx,nr,ny))
  DXygivenx <- array(0,c(nx,nr,ny))
  ## result of line below is p=nx by weight=nr by y=ny tensor.  Weight indexes rows of ygivenx, p indexes what deriv's are wrt.
  for (iw in 1:nr) {
    DPygivenx[ , iw, ] <- DProweight[ , iw] %o% (-(1/roweight[iw]^2) * eaUh[ipplus, ,drop=FALSE][iw, ])
  }
  for (iy in 1:ny) {
    DPygivenx[ , , iy] <- DPygivenx[ , , iy] + DPh[ , iy] %o% (c(1/roweight) * eaU[ipplus , iy])
  }
  for (iw in 1:nr) {
    DXygivenx[ , iw, ] <- DXroweight[ , iw] %o% ((-1/roweight[iw]^2) * eaUh[ixp[iw], ])
    DXygivenx[ixp[iw], iw, ] <- DXygivenx[ixp[iw],iw,] + (1/roweight[iw]) * eaUh[ixp[iw], ,drop=FALSE] * DXU[ixp[iw], ,drop=FALSE] *alph 
  }
  for (iy in 1:ny) {
    DXygivenx[ , , iy] <- DXygivenx[ , , iy] + DXh[ , iy] %o% (c(1/roweight) * eaU[ipplus , iy])
  }                                                           
  obj <- sum(f  * Umat ) - (1/alph) * pplus %*% sy(entel(ygivenx)) # + (1/alpha) * sum(gy * log(gy)) (doesn't change with param)
  DPobj <- c(DPobj) - c((1/alph) * (DPpnew[ , ipplus,drop=FALSE] %*% sy(entel(ygivenx)) + c(tensor(DPygivenx, c(pnew[ipplus]) * (1 + log(ygivenx)), 2:3, 1:2))))
  DXobj <- c(DXobj) - c((1/alph) * (DXpnew[ , ipplus,drop=FALSE] %*% sy(entel(ygivenx)) + c(tensor(DXygivenx, c(pnew[ipplus]) * (1 + log(ygivenx)), 2:3, 1:2))))
  obj <- -obj                           #as input to minimizer
  attr(obj, "gradient") <- -c(cbind(diag(nx - 1), rep(-1, nx - 1)) %*% DPobj, DXobj) #recognizing that there are only nx-1 p parameters, with  
  info = sum(f[f>0] * log(f[f>0])) - sum(log(pnew[pnew>0])*pnew[pnew>0]) - sum(log(gy[gy>0])*gy[gy>0])
  return(list(obj=obj, pnew=pnew, ygivenx=ygivenx, h=h, f=f, info=info, Eloss = sum(f * Umat)))
}
