infoD <- function(x,F,g,dh,dm){
  ## takes approximate model fitted with optDist and constructs diagnostics
  nm <- dim(F)[2]
  nh <- dim(F)[1]
  pm <- poly(1:nm,dm)
  pm <- cbind(rep(1/sqrt(nm),nm),pm)
  ph <- poly(1:nh,dh)
  ch <- x[1:dh]
  cm <- x[(dh+1):(dh+dm+1)]
  m <- exp(pm %*% cm)
  h <- exp(ph %*% ch)
  h <- h/sum(h)
  ones <- F %*% m
  ghat <- crossprod(h, F)*c(m)
  return(list(h=h,m=m,ones=ones,ghat=ghat))
}
