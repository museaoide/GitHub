g0g1dc <- function(ex,x,xd=paste(x,"dot",sep=""),shock){
  ## ex:     vector of expressions defining equilibrium
  ## x:      names (character) of current value variables 
  ## xd:     time-derivative  names. By default, just x names with "dot" appended.
  ## shock:  exogenous disturbance variable names
  nf <- length(ex)
  g0g1out <- vector("expression",length(ex))
  for(ix in 1:nf){
    g0g1out[[ix]] <- deriv(ex[ix],c(xd,x,shock))
  }
  names(g0g1out) <- names(ex)
  ## What is returned is a list of nf expressions.  When evaluated, these
  ## expressions yield the values of the expressions in ex, but they also have "gradient" attributes, so
  ## that attr(evaluated expression,"gradient") is the gradient vector,
  ## or a gradient matrix if the elements of x (say) are vectors.  Note that when this routine
  ## returns a matrix as the gradient, this is not the df/dx matrix with f the ex vector.  It is
  ## df(ix)/dx for a particular ix, evaluated at a vector of arguments.
  return(g0g1out)
}
