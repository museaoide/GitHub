numHess2 <- function(fcn, x, delta=NULL,...){
# delta: a scalar or a vector of increments. A scalar is applied to all the parameters
#        and increments of a vector are applied to corresponding parameters
	nx <- length(x)
	h <- matrix(0,nx,nx)
	f0 <- fcn(x,...)
	realsmall <- 0.00012207031250 # =.Machine$double.eps^(.25)
	if (is.null(delta)) { # rule of thumb
	  delta <- pmin(.5*abs(x),pmax(realsmall*abs(x),realsmall))
	}	
	nd <- length(delta)
	if ((nd==1)||(nd==nx)) {
	  tvec <- diag(delta,nrow=nx)
	} else {
		cat("Check the dimension of delta\n")
		return
	}
	for (indi in 1:nx){
		tvecv <- tvec[,indi]
    cat(sprintf("[%2.0f] [f+] %12.6f %5s [f-] %12.6f %5s [f0] %12.6f \n",indi,fcn(x+tvecv,...),(fcn(x+tvecv,...) > f0),fcn(x-tvecv,...),(fcn(x-tvecv,...) > f0),f0))
		h[indi,indi] <- (fcn(x+tvecv,...)+fcn(x-tvecv,...)-2*f0)/delta[indi]^2
	}
	if (nx > 1){
		for (indi in 1:(nx-1)){
			for (indj in (indi+1):nx){
				tvecv <- tvec[,indi] + tvec[,indj]
				h[indi,indj] <- .5*((fcn(x+tvecv,...)+fcn(x-tvecv,...)-2*f0)/delta[indi]^2 - h[indi,indi] - h[indj,indj])
				h[indj,indi] <- h[indi,indj]
			}
		}
	}
return(h)
}