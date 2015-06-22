mleBC <- function(ydata, lags, SigFacs, mnprior=list(tight=.2,decay=.5),vprior=list(sig=1,w=1), mu=10, Lambda, Mu, OmFac, cbar, csigFac=NULL) {
  require(Matrix)
  nT <- dim(ydata)[[1]]
  nv <- dim(ydata)[[2]]
  nc <- dim(ydata)[[3]]
  dataBlocks <- list(ActualData=NULL, BpriorDummies=NULL, y0ciDummies=NULL, ciprior=NULL)
  vp <- varprior(nv=nv,nx=0,lags=lags, mnprior=mnprior,vprior=vprior)
  ndobs <- dim(vp$ydum)[1] - (length(vp$pbreaks)+1) * lags
  ydum <- rbind(vp$ydum[vp$pbreaks, ], vp$ydum[dim(vp$ydum)[1], ])
  ## line above works because each dummy observation takes up exactly lags+1 positions in vp$ydum.  With training
  ## sample component (as in rfvar3 or mgnldnsty), would need more logic.
  xdum <- vp$ydum[-c(vp$pbreaks, dim(vp$ydum)[1]), ]
  xdum <- array(xdum, c(lags, ndobs, nv))
  xdum <- xdum[ seq(lags, 1, by=-1), , ]
  xdum <- aperm(xdum, c(2,1,3))
  xdum <- matrix(xdum, nrow=ndobs)
  ## Note that vp$xdum is exogenous variable data only.  xdum in this program is used to build whole
  ## rhs variable matrix, including lagged y's.  vp$xdum is not used at all in this program.
  ##------------------
  ## now sum of coeffs dummy obs
  ## Unlike rfvar3.R, here we don't use initial conditions for these.  This might also be better in rfvar3.
  ## Now the observations are scaled by vprior$sig, which means that the mu weight on them must be much bigger here
  ## than in rfvar3.R to have the same strength.  mu of 5 or 10 would be sensible.
  if (mu > 0) {
    ydumsc <- mu * diag(vprior$sig, nrow=nv)
    xdumsc <- mu * array(diag(vprior$si, nrow=nv), c(nv, nv, lags))
    xdumsc <- Matrix(as.vector(xdumsc), nrow=nv)
    ydum <- rbind(ydum, ydumsc)
    xdum <- rBind(Matrix(xdum), xdumsc)
  }
  ## get these prior dummies into regression matrix form
  dim(ydum) <- c(prod(dim(ydum)), 1)
  xdum <- kronecker(Matrix(1, nv, 1), xdum)  # stack up nv copies of xdum, one for each "variable" segment in ydum
  ## get actual data into regression matrix form
  y <- matrix(ydata[-(1:lags), , ], ncol=1)
  x <- array(0, c(nT-lags, nv, lags, nc))
  for (il in 1:lags)
    x[ , , il, ] <- ydata[(lags+1):nT - il, , ]
  X <- Matrix(0, nrow=0, ncol = nv * nv * lags)
  for (ic in 1:nc)
    X <- rBind(X, kronecker(Diagonal(nv), Matrix(x[ , , , ic], nrow = nT - lags))) #each variable in a country gets the same
                                        #matrix of rhs lagged y's
  ## now add all the country-specific constants as columns of X
  cdum <- array(diag(nv), c(nv, nv, nT-lags))
  cdum <- aperm(cdum, c(3,1,2))
  cdum <- Matrix(cdum, ncol=nv)
  X <- cBind(X, kronecker(Diagonal(nc), cdum))
  ## incorporate varprior and sum of coeff dummies into y and X
  y <- rbind(y, ydum)
  nActualRows <- dim(X)[1]
  dataBlocks$ActualData <- c(1,nActualRows)
  X <- rBind(X, Matrix(0, nrow=dim(xdum)[1], ncol=dim(X)[2]))
  X[-(1:nActualRows), 1:dim(xdum)[2]] <- xdum
  dataBlocks$BpriorDummies <- c(nActualRows+1, dim(X)[1])
  ##
  ## Now add rows that connect constants to y0's
  y0 <- Matrix(ydata[1:lags, , ], ncol=1) - Matrix(Mu, lags * nv * nc, ncol=1) #subtracts out the mu from y0=Lambda*cs+mu
  x0 <- Matrix(0, lags*nv*nc, dim(X)[2])
  x0[1:(lags * nv * nc) , lags * nv *nv + 1:(nc * nv)] <- kronecker(Diagonal(nc), Lambda)
  X <- rBind(X, x0)
  y <- rBind(y,y0)
  dataBlocks$y0ciDummies <- c(max(unlist(dataBlocks))+1,dim(X)[[1]])
  ## prior on country-variable constants
  if (!is.null(csigFac)) {
    y <- rBind(y, matrix(cbar, nv * nc, 1))
    xc <- Matrix(0, nv*nc, dim(X)[2])
    xc[ 1:(nv * nc), dim(X)[2] - nv * nc + 1:(nv * nc)] <- Diagonal(nv*nc)
    X <- rBind(X, xc)
    dataBlocks$ciprior <- c(max(unlist(dataBlocks))+1,dim(X)[[1]])    
  }
  ## Set up residual covariance matrix factor
  ## First actual data
  BigFac <- Matrix(0, length(y), length(y))
  for (ic in 1:nc) {
    icndx <- (ic - 1) * nv * (nT-lags) + 1:(nv * (nT-lags))
    BigFac[icndx, icndx] <- kronecker(SigFacs[ , , ic], Diagonal(nT-lags))
  }
  ## Then for prior dummy observations
  BpdNdx <- with(dataBlocks, BpriorDummies[1]:BpriorDummies[2])
  BigFac[BpdNdx, BpdNdx] <- kronecker(Diagonal(x=1/vprior$sig, n=nv), Diagonal(n=length(BpdNdx)/nv))
  ## Then y0 to ci relations
  y0ciNdx <- with(dataBlocks, y0ciDummies[1]:y0ciDummies[2])
  BigFac[y0ciNdx, y0ciNdx] <- kronecker(Diagonal(n=nc), OmFac)
  if (!is.null(csigFac)) {
    csigNdx <- with(dataBlocks, ciprior[1]:ciprior[2])
    BigFac[csigNdx, csigNdx] <- kronecker(Diagonal(nc), csigFac)
  }
  qrX <- qr(BigFac %*% X)
  stdy <- BigFac %*% y
  cf <- qr.coef(qrX, stdy)
  stdresid <- qr.resid(qrX, stdy)
  resid <- y - X %*% cf
  By <- array(cf[1:(nv * lags * nv)], c(nv, lags, nv)) #variables by lags by equations (order of X columns)
  By <- aperm(By, c(3,1,2))             #now equations by variables by lags, to match other VAR programs
  cs <- Matrix(cf[-(1:(nv * lags * nv))], nv, nc)
  cs <- t(cs)
  resid.std <- vector("list",0)
  resid.std$Actual <- matrix(as.vector(stdresid[dataBlocks$ActualData[1]:dataBlocks$ActualData[2]]), ncol=nv)
  resid.std$Bprior <- matrix(as.vector(stdresid[dataBlocks$BpriorDummies[1]:dataBlocks$BpriorDummies[2]]), ncol=nv)
  resid.std$y0ci <- matrix(as.vector(stdresid[dataBlocks$y0ciDummies[1]:dataBlocks$y0ciDummies[2]]), ncol=nv)
  resid.std$ciprior <- matrix(as.vector(stdresid[dataBlocks$ciprior[1]:dataBlocks$ciprior[2]]), ncol=nv)
  resid.raw <- vector("list",0)
  resid.raw$Actual <- matrix(as.vector(resid[dataBlocks$ActualData[1]:dataBlocks$ActualData[2]]), ncol=nv)
  resid.raw$Bprior <- matrix(as.vector(resid[dataBlocks$BpriorDummies[1]:dataBlocks$BpriorDummies[2]]), ncol=nv)
  resid.raw$y0ci <- matrix(as.vector(resid[dataBlocks$y0ciDummies[1]:dataBlocks$y0ciDummies[2]]), ncol=nv)
  resid.raw$ciprior <- matrix(as.vector(resid[dataBlocks$ciprior[1]:dataBlocks$ciprior[2]]), ncol=nv)
  return(list(By=By, cs=cs, resid.std=resid.std, resid.raw=resid.raw, qrX=qrX, stdy=stdy, dataBlocks=dataBlocks) )
}
