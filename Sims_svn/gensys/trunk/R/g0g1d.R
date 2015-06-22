#' g0g1d
#' 
#' Find analytic derivatives of an \code{eqsys} object.
#' 
#' @param ex \code{eqsys} object --- an equation system.
#' @param x  character vector of names of current value variables 
#' @param xl character vector of lagged value names. By default, just \code{x}
#'   names with "l" appended.
#' @param shock  character vector of exogenous disturbance variable names
#' @return
#' A list of nf expressions.  When evaluated, these
#' expressions yield the values of the expressions in ex, but they also
#' have "gradient" attributes, so that
#' \code{attr(evaluated expression,"gradient")} is, for each of the
#' expressions in \code{ex},  the gradient vector for that \code{ex[[ix]]}.
#'
#' @details
#' This function produces expressions, not values.  To get numeric values for
#' the derivative matrices, the output from this function is handed to
#' \code{g0g1eval}, along with the steady state value (from \code{ssSolve()})
#' and parameter settings to use in doing the evaluation.
#'
#' @export
g0g1d <- function(ex, x=attr(ex,"vlist"), xl=paste(x,"l",sep=""), shock=attr(ex, "shock")) {
  nf <- length(ex)
  if(shock[1] == "NONE") shock <- vector("character",0)
  g0g1out <- vector("expression",length(ex))
  for(ix in 1:nf){
    g0g1out[[ix]] <- deriv(ex[ix],c(x,xl,shock))
  }
  g0g1out <- structure(g0g1out, class="eqsys", names=names(ex), vlist=attr(ex, "vlist"), shock=attr(ex, "shock"), param=attr(ex, "param"), forward=attr(ex,"forward"))
  return(g0g1out)
}
