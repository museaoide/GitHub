#' ssSolve
#'
#' Find the steady state of a dynamic \code{eqsys} object.
#'
#' @param ex the \code{eqsys} object
#' @param x0 an initial guess of the steady state
#' @param crit size of the discrepancy in steady state equations that
#'   is treated as effectively zero.
#' @param itmax maximum number of iterations of the nonlinear solver.  The
#'   default 20 is usually enough, but up to 400 is also usually fast.
#' @param Should the iterations of the solver be printed out in full?
#' @param alpha,delta See the documentation for \code{csolve()}
#' @param long If \code{verbose}, do not include full printout of
#'   function and parameter value vectors.
#'
#' @return
#' \describe{
#'    \item{xss} steady state value of the variables
#'    \item{csout} returned object from the nonlinear equation solver
#'    \item{fss} the function, constructed from ex, that returns zero
#'       when its argument is a steady state value
#' }
#' @export
ssSolve <- function(ex, x0, param, crit=1e-7, itmax=20, verbose=TRUE, alpha=1e-3, delta=1e-6, long=FALSE) {
  ex2 <- sssys(ex)
  shocknames <- attr(ex,"shock")
  shockval <- rep(0, length(shocknames))
  names(shockval) <- shocknames
  fss <- derivVec(ex2, attr(ex2, "vlist"), c(param, shockval))
  ## fss returns the value of the ex2 system expressions, and also their derivatives.  csolve() knows how to use such a function.
  csout <- csolve(fss, x0, gradfun=fss, crit=crit, itmax=itmax, verbose=verbose, alpha=alpha, delta=delta, long=long)
  xss <- structure(as.vector(csout$x), names=dimnames(csout$x)[[1]], param=param)
  return(list(xss=xss, csout=csout, fss=fss))
}
