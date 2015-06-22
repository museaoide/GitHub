ssSolve <- function(ex, x0, param, crit=1e-7, itmax=20, verbose=TRUE, alpha=1e-3, delta=1e-6, long=FALSE) {
  ex2 <- sssys(ex)
  shocknames <- attr(ex,"shock")
  shockval <- rep(0, length(shocknames))
  names(shockval) <- shocknames
  fss <- derivVec(ex2, attr(ex2, "vlist"), c(param, shockval))
  csout <- csolve(fss, x0, gradfun=fss, crit=crit, itmax=itmax, verbose=verbose, alpha=alpha, delta=delta, long=long)
  xss <- structure(as.vector(csout$x), names=dimnames(csout$x)[[1]], param=param)
  return(list(xss=xss, csout=csout, fss=fss))
}
