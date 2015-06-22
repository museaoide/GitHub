ssf <- function(xval,dfcn,xdotval,shock,experr,param){
  ## xdotval: named list of values for xdot at which we are evaluating the function and gradient. (all 0)
  ## xval:    named list of values for x at which we are evaluating the function and gradient. 
  ## shock:  named list of exogenous shocks (all 0)
  ## experr:   named list of expectational errors (all 0)
  ## param:   named list of parameter values at which we are evaluating the function and gradient
  ## dfcn:    the returned list from g0g1
  ##-------------------
  ## discrep: the vector of values of the expression list that entered g0g1 as ex
  ## grad:    the nv x nf (square) matrix of derivatives of the discrepancy w.r.t. the x vector
  ##-------------------
  ## The "all 0" notes apply to use of this as input to solving for steady-state.  One could in principle
  ## use it to solve for xval with particular fixed shock values and fixed xdotval that are not zero.
  nv <- length(xval)
  nf=length(dfcn$g0)
  ## if(!(nv==nf))stop("mismatch of x and f lengths") # wrong place to check this
  discrep <- matrix(0,nf,1)
  grad <- matrix(0,nf,nv)
  for(ifcn in 1:nf){
    f <- dfcn$g1[[ifcn]]
    dball <- do.call("f",c(xdotval,xval,shock,experr,param))
    discrep[ifcn,1] <- c(dball)
    grad[ifcn,] <- attr(dball,"gradient") # for some reason the column names in gradient are lost here
  }
  dimnames(discrep) <- list(names(dfcn$g1),list(NULL))
  browser()
  dimnames(grad) <- list(names(dfcn$g1),names(xval))
  return(list(discrep=discrep,grad=grad))
}
