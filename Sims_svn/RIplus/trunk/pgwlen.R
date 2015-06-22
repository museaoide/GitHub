pgwlen <- function(nx,nw,xrange,wrange) {
  x <- seq(xrange[1],xrange[2],length.out=nx)
  w <- seq(wrange[1],wrange[2],length.out=ny)
  cltw <- outer(x,w,function(a,b){a<b})
  return(sum(cltw[2:nx, ]))
}
