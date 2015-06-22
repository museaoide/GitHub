screpNew <- function(params=vector(numeric),nmix=3,U=function(a,b){log(ifelse(a <= 0,1e-100,a))+log(ifelse(b-a <= 0,1e-100,b-a))},
                     xrange=c(0,4),nx,gg=matrix(0,800,2),lambda=1)
{
  ## parameterizes h(x) as a mixture of nmix normals with means mu and common variance sigma.
  ## Tried also mixture of gammas with equal scale
  ## and most recently straight spline, with negative elements set to zero
  ## Scales h(x)*exp(lambda*U(x,w))*gg(w) to integrate to one.
  ## Forms: EU = sum over (x,w) of h(x)*exp(lambda*U(x,w))gg(w)*U(x,w) (expected utility)
  ##        scvec = sum over x of h(x)exp(lambda*U(x,w))-1, for each w (vector of discrepancies with FOC)
  ##        scprj = vector of crossprods of scvec with polynomials
  ##        obf = EU - sum(scvec^2)*scweight (penalized objective function
  ## params: vector packed with parameters (mu, sigma, weights)
  ## nmix:  number of shifted normals used
  ## U:     function(x,w), the utility function
  ## x:     equispaced grid of x values
  ## gg:     the given marginal pdf of w; an n x 2 vector, with first column the pdf values, 2nd column the (equispaced) w values.  sum(gg[,1]==1.
  ## scweight:  Weight on the constraint in the penalized objective function, in non-spline versions
  ## lambda:inverse of Lagrange multiplier on info constraint.  Increases with increasing capacity.
  npar <- length(params)
  ## mu <- params[1:nmix]
  ## mu <- exp(params[1:nmix])
  ## mu[mu<1e-10] <- 1e-10
  ## mu[mu>1e10] <- 1e10
  ## sig <- exp(params[nmix+1])
  ## w <- params[(nmix+2):(2*nmix)]
  ##-------------
  ## hx is 0 and 1 plus points in between, with in-between points the first (npar-1)/2 parameters.
  ## hy are the postulated values of hy at the hx points, normalized to sum to one.
  hx <- c(0,params[1:((npar-1)/2)],1)*(xrange[2]-xrange[1])+xrange[1]
  hy <- params[((npar+1)/2):npar]
  hy <- hy/(1+sum(hy))
  hy <- c(hy,1-sum(hy))
  h <- spline(hx,hy,n=nx)
  h$y <- pmax(h$y,0)
  h$y <- h$y/sum(h$y)
  ##-----------------
  ## x <- as.array(x)
  ## scprj <- matrix(0,npar,1)
  ## x <- seq(-.5,.5,.01)
  ##if(any(w>200)){
  ##  d1 <- (w[1]>200)
  ##  d2 <- w[2]>200
  ##  d1 <- d1/(d1+d2)
  ##  d2 <- 1-d1
  ##  d3 <- 0
  ##}else{
  ##  dt <- 1+exp(w[1])+exp(w[2])
  ##  d1 <- exp(w[1])/dt
  ##  d2 <- exp(w[2])/dt                  
  ##  d3 <- 1-d1-d2
  ##}
  W <- gg[,2,drop=FALSE]
  ## xs <-drop(outer(as.array(x),as.array(1:npar),"^"))
  ## xs <- qr.Q(qr(xs))
  ## h <- as.array(d1*dnorm(x,mu[1],sig)+d2*dnorm(x,mu[2],sig)+d3*dnorm(x,mu[3],sig))
  ## h <- as.array(d1*dgamma(x,shape=mu[1],scale=sig)+d2*dgamma(x,shape=mu[2],scale=sig)+d3*dgamma(x,shape=mu[3],scale=sig))
  ## h <- h/sum(h)
  Umat <- drop(outer(h$x,W,FUN=U))
  pdf <- h$y * exp((Umat-max(Umat))*lambda)
  pdf[drop(outer(h$x,W,">"))] <- 0
  m <- gg[,1]/apply(pdf,MAR=2,FUN=sum)
  pdf <- t(m * t(pdf))
  pdf[!is.finite(pdf)] <- 0
  EU <- sum(pdf * Umat)
  ##  scvec <-  exp(Umat) %*% as.vector(m) - 1 
  ##  if(any(is.nan(scvec)) || max(abs(scvec)) > 1e100){
  ##    scprj[] <- 1e200
  ##    EU <- -1e200
  ##    obf <- -1e200
  ##    scvec[] <- 1e200
  ##  }else{
  xpdf <- apply(pdf,1,sum)
  obf <- EU-(1/lambda)*(sum(log(pdf[pdf>0])*pdf[pdf>0])-sum(log(xpdf[xpdf>0])*xpdf[xpdf>0]))
  ## obf <- EU - scweight*sum(scvec^2)
  ## scvec is now length(x).  xs is length(x) x npar(indexing powers)
  ##  scprj <- apply(as.vector(scvec) * xs, MAR=2,FUN=sum)
  ##  scprj <- t(scprj)
  ## return(list(obf=obf,EU=EU,scprj=scprj,scvec=scvec,pdf=pdf))
  return(list(obf=obf,EU=EU,pdf=pdf,h=h))
##return(-obf) # For use in csminwel
}

screpCS <- function(...){-screpNew(...)$obf}
