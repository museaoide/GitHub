makepdf <- function(xh,nx,nw,xrange,wrange) {
  x <- seq(xrange[1],xrange[2],length.out=nx)
  w <- seq(wrange[1],wrange[2],length.out=nw)
  ggw <- gg(w)
  ggw <- ggw/sum(ggw)
  cltw <- outer(x,w,function(a,b){a<b})
  pgwLoc <- matrix(-100,nx-1,nw)
  pgwLoc[cltw[2:nx,]] <- xh
  ##-------------------------------------------------------------------------------
  ## Map pgw into a matrix whose values are all in (0,1) and sum to one in each col.
  pgwLoc <- exp(pgwLoc)
  pgwLoc <- rbind(rep(1,nw), pgwLoc)
  colsumi <- 1/(as.vector(rep(1,nx) %*% pgwLoc))
  pgwLoc <- t(colsumi * t(pgwLoc))
  ##--------------------------------------------------------------------------------
  return(pdf <- t(ggw * t(pgwLoc)))
}
