makessfct <- function(eq, xnames, dxnames=paste(xnames,"dot",sep=""), shocknames, paramnames, paramv) {
  ## creates a steady-state function with a single numeric matrix argument, with the fixed parameter
  ## and (zero) shock values stored in an environment within this makessf function.
  eqev <- vector("numeric", 2*length(xnames)+length(shocknames)+length(paramnames))
  names(eqev) <- c(xnames,dxnames,shocknames,paramnames)
  eqev[paramnames] <- paramv
  eqev[shocknames] <- 0
  ssf <- function(xv) {
    xv <- as.matrix(xv)
    nx <- dim(xv)[2]
    fval <- matrix(0,length(eq),nx)
    for (ix in 1:nx) {
      eqev[xnames] <- xv[ , ix]
      eqev[dxnames] <- 0
      for (iq in 1:length(eq)) {
        fval[iq, ix] <- eval(eq[iq],as.list(eqev))
      }
    }
    return(fval)
  }
  return(ssf)
}
