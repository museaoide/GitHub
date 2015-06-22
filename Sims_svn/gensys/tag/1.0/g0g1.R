g0g1 <- function(ex,x,xl,shock,experr,param){
  ## ex:     vector of expressions defining equilibrium
  ## x:      current values (or time derivatives in continuous case)
  ## xl:     lagged values (or levels in continuous case)
  ## shock: exogenous disturbances
  ## experr: expectational errors
  ## param:  other arguments, with respect to which we are not differentiating
  nf <- length(ex)
  g0 <- list(0)
  g1 <- list(0)
  Phi <- list(0)
  Pi <- list(0)
  for(ix in 1:nf){
    g0[[ix]] <- deriv(ex[ix],x) #,function.arg=c(x,xl,shock,experr,param) )
    g1[[ix]] <- deriv(ex[ix],xl) #,function.arg=c(x,xl,shock,experr,param))
    Phi[[ix]] <- deriv(ex[ix],shock) #,function.arg=c(x,xl,shock,experr,param))
    Pi[[ix]] <- deriv(ex[ix],experr) #,function.arg=c(x,xl,shock,experr,param))
  }
  names(g0) <- names(ex)
  names(g1) <- names(ex)
  names(Phi) <- names(ex)
  names(Pi) <- names(ex)
  ## What is returned is a list of 4 lists of nf expressions.  When evaluated, these
  ## expressions yield the values of the expressions in ex, but they also have "gradient" attributes, so
  ## that attr(evaluated expression,"gradient") is the gradient vector,
  ## or a gradient matrix if the elements of x (say) are vectors.  Note that when this routine
  ## returns a matrix as the gradient, this is not the df/dx matrix with f the ex vector.  It is
  ## df(ix)/dx for a particular ix, evaluated at a vector of arguments.
  return(list(g0=g0,g1=g1,Phi=Phi,Pi=Pi))
}
