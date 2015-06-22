makessf <- function(eq, xnames, xlnames=paste(xnames,"l",sep=""), shocknames, paramnames, paramv){
  ## creates a steady-state function with a single numeric vector argument, with the fixed parameter
  ## and (zero) shock values stored in an environment within this makessf function.
  eqev <- vector("numeric", 2*length(xnames)+length(shocknames)+length(paramnames))
  names(eqev) <- c(xnames,xlnames,shocknames,paramnames)
  eqev[paramnames] <- paramv
  eqev[shocknames] <- 0
  eqev <- as.list(eqev)
  ssf <- function(xv) {
    fval <- vector("numeric",length(eq))
    eqev[xnames] <- xv
    eqev[xlnames] <- xv
    for (iq in 1:length(eq)) {
      fval[iq] <- eval(eq[iq],eqev)
    }
    return(fval)
  }
  return(ssf)
}
## this seems to be a "pre-eqsys" version of steady state code.
