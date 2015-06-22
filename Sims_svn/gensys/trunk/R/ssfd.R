ssfd


#' xval:    named list of values for x at which we are evaluating the function and gradient. 
#' shock:  named list of exogenous shocks (all 0)
#' param:   named list of parameter values at which we are evaluating the function and gradient
#' dfcn:    the returned list from g0g1d
#' xlval:  named list of values for xl at which we are evaluating the function and gradient. Default
#'         just xval, with "l" added to all the names.
#'-------------------
#' discrep: the vector of values of the expression list that entered g0g1d as ex
#' grad:    the nv x nf (square) matrix of derivatives of the discrepancy w.r.t. the x vector
#'-------------------
ssfd <- function(xval, dfcn, shock, param, xlval=xval){
  if(identical(xval,xlval)) {
    names(xlval) <- paste(names(xval),"l",sep="l")
  }
  nv <- length(xval)
  nf=length(dfcn$g0)
  ## if(!(nv==nf))stop("mismatch of x and f lengths") # wrong place to check this
  discrep <- matrix(0,nf,1)
  grad <- matrix(0,nf,nv)
  
  dimnames(discrep) <- list(names(dfcn$g1),list(NULL))
  dimnames(grad) <- list(names(dfcn$g1),names(xval))
  return(list(discrep=discrep,grad=grad))
}
