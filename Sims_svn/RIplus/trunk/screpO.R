screpO <- function(params){
  ## total discrepancy, for minimization.
  npar <- length(params)
  nn <- npar/2
  mu <- params[1:nn]
  sig <- exp(params[nn+1])
  w <- params[(nn+2):npar]
  x <- seq(-.5,.5,.01)
  nx <- length(x)
  topw <- max(w)
  if(topw>100) {                        # no weight on first
    w <- w-topw
    dt <- sum(exp(w))
  } else {
    dt <- 1+sum(exp(w))
  }
  d <- exp(w)/dt
  d <- c(d,1-sum(d))
  x <- rep(x,each=nn)
  pdf <- as.vector(d)*dnorm(x,mu,sig)  # mu is recycled by dnorm, while d is recycled by the product
  dim(pdf) <- c(nn,nx)
  ## dim(lpdf) <- c(nn,nx)
  ## lpdf <- log(apply(lpdf,2,sum))
  ## return(sum(lpdf^2))
  pdf <- apply(pdf,2,sum)
  return(-sum(log(pdf)))
}

            
