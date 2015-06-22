screpP <- function(params){
  x <- seq(-1,1,.01)
  np <- length(params)
  if(params[length(params)]<0){
    return(1e100)
  }
  fh <- exp(-( outer(x,2*(1:np),"^") %*% params ))
  fh <- 101*fh/sum(fh)
  return(sum(abs(fh[51:151]-1))+sum(abs(fh[1:50]))+sum(abs(fh[152:201])))
}
