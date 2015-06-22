"gstate" <-
function(G1,impact,pickS=NULL,pickC=NULL) {
### G1:    the coefficient on  lagged y from the output of gensys.m.  
### impact:the coefficient on exogenous shocks from the output of gensys.m
### pick:  an optional guess at a matrix of coefficients that extracts 
###        the state vector from y.  Must be a matrix with nrow < ncol
### ok:     0 if pick matrix is just no good.
###         1 if pick matrix is usable forward, but loses initial date information
###         2 if pick matrix is not usable forward, is part of a correct state vector, but is not complete
###         3 if pick matrix is     usable forward, is part of a correct state vector, but is not complete
###         4 if pick matrix is not usable forward, summarizes past, but is redundant
###         5 if pick matrix is     usable forward, summarizes past, but is redundant
###         6 if pick matrix summarizes past and is not redundant, but is not usable forward
###         7 if pick matrix is perfect, both forward and as history summary
### pickSn: a matrix of coefficients that extracts the state vector from y.  Equal
###        to pick if pick is supplied and ok is odd; otherwise a matrix similar to pick that generates an odd ok.
### GS:    the matrix of coefficients on the lagged state variable.
### uis,vs: If uis'*vs is full rank, ok=7 is possible, otherwise not.  If ok<2, an ok>2 can be found by trying
###        a pick in the row space of vs'.  Any pick with pick*uis full column rank will provide a foward state 
###        (i.e. ok an odd number).  {<- seems not}
### The solution was in the form y(t)=G1 %*% y(t-1)+impact %*% z(t).  Now it's in the form y(t)=GS %*% pickn %*% y(t-1)+impact %*% z(t).
### In general, pickn is mxn with m<n.
  REALSMALL <- 1e-9
  nr <- dim(G1)[1]
  if(dim(G1)[1]!=dim(G1)[2]) stop("G1 must be square")
  ##if nr<nc
  ##   G1=[G1;zeros(nc-nr,nc)]
  ##end
  svdG1 <- svd(G1)
  top <- svdG1$d > REALSMALL
  nd <- sum(top)
  us <- svdG1$u[,top,drop=FALSE]
  uu <- svdG1$u[,!top,drop=FALSE]
  vs <- svdG1$v[,top,drop=FALSE]
  if (is.null(pick)){
    pick <- t(vs)
    vp <- svdG1$v
    vps <- vs
    dp <- rep(1,nd)
    ups <- diag(nrow=nd,ncol=nd)
    ndp <- nd
  }else{
    svdp <- svd(pick)
    topp <- svdp$d > REALSMALL
    ndp <- sum(topp)
    vps <- svdp$v[,topp,drop=FALSE]
    dp <- svdp$d[topp]
    ups <- svdp$u[topp,topp,drop=FALSE]
  }
  ## If we were worried about efficiency, we'd skip some of this when pick=vs.
  ##Does pick summarize history?  (pick in v' row space, v in vp' row space).
  pinv <- all( (pick %*% (diag(nr)-vs%*%t(vs)))^2 < REALSMALL ) 
  vinp <- all(((diag(nr)-vps%*% t(vps)) %*% vs)^2 < REALSMALL)
  okpast <- pinv && vinp
  ## Does pick summarize all current info?  (pick spans column space of cbind(impact,G1) )
  svdi <- svd(cbind(impact, us))
  topi <- svdi$d > REALSMALL
  ndi <- sum(topi)
  uis <- svdi$u[,topi,drop=FALSE]
  uiu <- svdi$u[,!topi,drop=FALSE]
  if (ndi < dim(G1)[1]) {
    if (dim(pick)[1] < dim(uis)[2]) {
      oknow <- FALSE
    } else {
      svdt <- svd(pick %*% uis)
      toppu <- svdt$d > REALSMALL        
      oknow <- sum(toppu) >= dim(uis)[2]
    }
  }else{
    oknow <- FALSE
  }
  if (vinp) {
    #browser()
    GS <- G1 %*% vps %*% ((1/dp)*t(ups))  #G1/pick 
    pickn <- pick
  } else {
    if (pinv) {
      r <- vs-vps %*% t(vps) %*% vs
      svdr <- svd(r)
      topr <- svdr$d > REALSMALL
      p2 <- svdr$u[,topr,drop=FALSE]
      pickn <- rbind(pick,t(p2))
      svdpn <- svd(pickn)
      GS <- G1 %*% svdpn$v %*% ((1/svdpn$d)*t(svdpn$u))  # G1/pickn
    } else {
      if (oknow) {
        GS <- G1 %*% solve(rbind(pick,t(uiu)))
        GS <- GS[,1:dim(pick)[1]]
        pickn <- pick
      } else {
        pickn <- t(vs)
        GS <- t(svdG1$d[top] * t(us))
      }
    }
  }
  ok <- oknow+2*pinv+4*vinp
  return(list(GS=GS,pickn=pickn,ok=ok,uis=uis,vs=vs))
}

