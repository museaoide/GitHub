g0g1d <- function(ex, x=attr(ex,"vlist"), xl=paste(x,"l",sep=""), shock=attr(ex, "shock")) {
  ## ex:     vector of expressions defining equilibrium
  ## x:      names (character) of current value variables 
  ## xl:     lagged value names. By default, just x names with "l" appended.
  ## shock:  exogenous disturbance variable names
  nf <- length(ex)
  g0g1out <- vector("expression",length(ex))
  for(ix in 1:nf){
    g0g1out[[ix]] <- deriv(ex[ix],c(x,xl,shock))
  }
  names(g0g1out) <- names(ex)
  g0g1out <- structure(g0g1out, names=names(ex), vlist=attr(ex, "vlist"), shock=attr(ex, "shock"), param=attr(ex, "param"), forward=attr(ex,"forward"))
  ## What is returned is a list of nf expressions.  When evaluated, these
  ## expressions yield the values of the expressions in ex, but they also have "gradient" attributes, so
  ## that attr(evaluated expression,"gradient") is, for each of the expressions in ex,  the gradient vector for that ex[[ix]].
  return(g0g1out)
}
