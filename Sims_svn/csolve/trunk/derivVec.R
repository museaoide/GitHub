#' derivVec
#'
#' Application of \code{deriv()} to a vector
#' 
#' @param ex  a vector of expressions, or an \code{eqsys} object
#' @param x    a character vector of variable names w.r.t. which the expressions
#'             will be differentiated.  If \code{ex} is an \code{eqsys} object,
#'             the default is for \code{x} to be the \code{vlist} attribute of
#'             \code{ex}.  
#' @param param  a vector of named parameter values that will be constant in
#'               repeated calls to \code{fret} (the returned function)
#' @param xchk    Range check function.  It returns \code{FALSE} if given an
#'                \code{x} vector that is outside the domain of definition of
#'                \code{ex}.
#'-------------------
#' @return \code{fret}, a function that when evaluated at a numerical \code{x}, 
#'           returns the vector of expression values, but also, as
#'           \code{attr(value, "grad")}, the gradient matrix
#'           It uses the fixed \code{param} values set in the call to
#'           \code{derivVec}.
#' @export
#' 
derivVec <- function(ex, x=attr(ex,"vlist"), param=vector("numeric",0), xchk=function(z){TRUE}) {
  nq <- length(ex)
  nv <- length(x)
  param <- param #to put param into this namespace, so it does not need to be in every call.
  xchk <- xchk
  outf <- vector("expression",nq)
  for (iq in 1:nq) {
    f <- deriv(ex[iq], x)
    outf[iq] <- f
  }
  fret <- function(z) {
    fval <- vector("numeric", nq)
    gval <- matrix(0, nq, nv)
      zv <- as.vector(z)
      names(zv) <- if (is.null(names(z))) dimnames(z)[[1]] else names(z)
    if (xchk(zv) ) {
      for (iq in 1:nq) {
        fg <- eval(outf[iq], as.list(c(zv, param)))
        fval[iq] <- fg
        gval[iq, ] <- attr(fg, "grad")
      }
      fval <- c(fval)
    } else {                            #bad z input
      fval[] <- 1e20
      gval <- diag(-1, nq, nv)
    }
    names(fval) <- names(ex)
    dimnames(gval) <- list(eq=names(ex), vbl=attr(ex,"vlist"))
    attr(fval, "grad") <- gval
    return(fval)
  }
  return(fret)
}
