infoDwel <- function(x,F,g){
  ## takes approximate model fitted with optDist and constructs diagnostics
  nm <- dim(F)[2]
  nh <- dim(F)[1]
  dh <- length(x)
  ph <- poly(1:nh,dh)
  h <- exp(ph %*% x)
  h <- as.vector(h/sum(h))
  m <- as.vector(as.vector(g)/(h %*% F))
  ones <- F %*% m
  hhat <- h * ones
  umat <- matrix(0,nh,nm)
  umat[F>0] <- log(F[F>0])
  U <- h %*% (F * umat) %*% m
  imat <- matrix(0,nh,nm)
  imat <- h * t(m * t(F))
  imat[imat>0] <- log(imat[imat>0])
  minfo <- h %*% (F * imat) %*% m -sum(g[g>0]*log(g[g>0])) - sum(hhat[hhat>0]*log(hhat[hhat>0]))
  ## browser()
  return(list(h=h,m=m,ones=ones,U=U,minfo=minfo))
}
