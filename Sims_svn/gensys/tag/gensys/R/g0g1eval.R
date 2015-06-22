#' g0g1eval
#'
#' Evaluates linearized system matrices at given parameters.
#' 
#'   \code{x}, \code{xl}, and \code{param} must all be vectors with named
#'   components, or lists.  When \code{x} is something like \code{flmss$xss},
#'   where \code{flmss} is a returned value from \code{ssSolve}, the default
#'   assignment for \code{param} works.
#'
#' @param dexpr   derivative expression vector, as produced by \code{g0g1d}
#' @param    x     vector of named values (not names) of \code{x}
#' @param   xl     vector of named values of \code{xl} 
#' @param shock    a vector of named values of exogenous shocks (all assumed
#'     centered at 0)
#' @param experr logical vector, TRUE for equations with expectational errors
#' @param param  vector of values of parameters (all other symbols) in
#'     the expressions.
#' 
#'
#' @return
#'   A list for use as inputs to \code{gensys()}, with elements
#'   \describe{
#'   \item{\code{g0}}{left-hand side coefficient matrix}
#'   \item{\code{g1}}{right-hand side coefficient matrix}
#'   \item{\code{Psi}}{matrix of coefficients on exogenous shocks}
#'   \item{\code{Pi}}{matrix of coefficients on endogenous expectation errors}
#'   \item{\code{param}}{parameter vector (same as input \code{param})}
#'   }
#'
#' @export
#' 
g0g1eval <- function(dexpr, x, xl=as.vector(x), shock=attr(dexpr, "shock"), experr=attr(dexpr, "forward"), param=attr(x,"param")) {
    ##             
    if (is.null(names(x))) { names(x) <- dimnames(x)[[1]]}
    if (is.null(names(xl))) names(xl) <- paste(names(x),"l",sep="")
    if(identical(shock, "NONE")) {
        shock <- vector("character", 0)
    } else {
        if(is.character(shock)) {
            nshock <- shock
            shock <- rep(0, length(nshock))
            names(shock) <- nshock
        } else {
              if (is.null(names(shock))) names(shock) <- dimnames(shock)[[1]]
          }
    }
    if (is.null(names(param))) names(param) <- dimnames(param)[[1]]
    nex <- length(dexpr)
    nx <- length(x)
    nxl <- length(xl)
    nshock <- length(shock)
    g0 <- matrix(0, nex, nx)
    g1 <- matrix(0, nex, nxl)
    Psi <- matrix(0, nex, nshock)
    nxerr <- sum(experr)
    Pi <- matrix(0, nex, sum(experr))
    for ( i in 1:nex ) {
        dvec <- attr(eval(eval(dexpr[i]), as.list(c(x, xl, shock, param))),"gradient")
        loc <- 0
        g0[i,] <- dvec[1:nx]
        loc <- nx
        g1[i,] <- dvec[(loc+1):(loc+nxl)]
        loc <- loc + nxl
        Psi[i,] <- dvec[(loc+1):(loc+nshock)]
    }
    Pi[experr,] <- diag(nxerr)
    dimnames(g0) <- list(names(dexpr),names(x))
    dimnames(g1) <- dimnames(g0)
    ## if(is.numeric(experr)) experr <- names(dexpr)[experr]
    dimnames(Psi) <- list(names(dexpr),names(shock))
    dimnames(Pi) <- list(names(dexpr), names(dexpr)[experr])
    return(list(g0=g0, g1=-g1, Psi=-Psi, Pi=-Pi, param=param))
}
