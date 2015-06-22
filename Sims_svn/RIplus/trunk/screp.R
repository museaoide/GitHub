screp <- function(params){
  ncol <- dim(params)[2]
  npar <- dim(params)[1]
  mu <- params[1:3,,drop=FALSE]
  sig <- exp(params[4,,drop=FALSE])
  w <- params[5:6,,drop=FALSE]
  s <- matrix(0,npar,ncol)
  x <- seq(-.5,.5,.01)
  nx <- length(x)
  if(any(w>200)){
    d1 <- (w[1,]>200)
    d2 <- w[2,]>200
    d1 <- d1/(d1+d2)
    d2 <- 1-d1
    d3 <- 0
  }else{
    dt <- 1+exp(w[1,])+exp(w[2,])
    d1 <- exp(w[1,])/dt
    d2 <- exp(w[2,])/dt
    d3 <- 1-d1-d2
  }
  xs <- outer(x,1:npar,"^")
  xs <- array(rep(xs,ncol),dim=c(nx,npar,ncol))
  xs <- aperm(xs,c(3,1,2))
  x <- matrix(rep(x,ncol),ncol=ncol)
  x <- t(x)
  lpdf <- log(d1*dnorm(x,mu[1,],sig)+d2*dnorm(x,mu[2,],sig)+d3*dnorm(x,mu[3,],sig))
  ## lpdf is now ncol x length(x).  xs is ncol x length(x) x npar(indexing powers)
  if(any(is.nan(lpdf)) || min(lpdf)< -1e10){
    s <- s - 1e20
  }else{
    s <- apply(xs*as.vector(lpdf),MAR=c(1,3),FUN=sum)
    s <- t(s)
  }
  ## browser()
  return(s)
}

            
