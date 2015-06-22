blmWrap <- function(par, x, y, dfp=1, scalepv, vlist, dvname="y", lags=rep(0, length(vlist)), ldv=1, scale, 
bbar=NULL, vmeans) {
  smooth <- par[1]
  damp <- par[2]
  erratio <- par[3]
  tsp <- tsregPrior(vlist, dvname, lags, ldv, scale, bbar,  smooth, damp, vmeans, erratio)
  blmout <- blm(x, y, tsp$x, tsp$y, dfp, scalepv)
  return(blmout$lmdd)
}