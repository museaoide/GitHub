dlogitbeta <- function(x,alpha,beta,location=0,scale=1,log=FALSE){
  q <- alpha/(alpha+beta)
  p <- 1-q
  z <- (x - location)/scale
  dn <- p^(-q)*q^(-p)/(exp(alpha*z)/alpha + exp(-beta*z)/beta)/scale/beta(a=q,b=p)
  if (log) return(log(dn)) else return(dn)
}

rlogitbeta <- function(n,alpha,beta,location=0,scale=1){
  q <- alpha/(alpha+beta)
  p <- 1-q
  z <- rbeta(n=n,shape1=p,shape2=q,ncp=0)
  x <- (log(alpha/beta) + log(z/(1-z)))/(alpha+beta)
  x <- scale*x + location
  return(x)
}

plogitbeta <- function(qntile,alpha,beta,location=0,scale=1,log.p=FALSE){
  q <- alpha/(alpha+beta)
  p <- 1-q
  xqntile <- (qntile - location)/scale
  zqntile <- beta/alpha*exp((alpha+beta)*xqntile) / (1 + beta/alpha*exp((alpha+beta)*xqntile))
  pr <- pbeta(q=zqntile, shape1=p, shape2=q, ncp = 0, lower.tail = TRUE, log.p=log.p)
  return(pr)
}