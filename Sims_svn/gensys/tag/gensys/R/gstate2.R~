gstate <-
function(G1, impact, pickS=NULL,pickC=NULL, allowRedundant=FALSE) {
### Inputs:
### G1:    the coefficient on  lagged y from the output of gensys.m.  
### impact:the coefficient on exogenous shocks from the output of gensys.m
### pickS:  an optional guess at a matrix of coefficients that extracts 
###        the state vector from y.  Must be a matrix with nrow < ncol
### pickC:  an optional guess at a matrix of coefficients that extracts 
###        the "control" vector from y.  Must be a matrix with nrow < ncol
### outputs:
### pickS: a matrix of coefficients that extracts a state vector from y.  Equal
###        to input pickS if okS is TRUE;
### pickC: a matrix of coefficients that extracts a "control" vector from y.  Equal to
###        input pickC if okC is TRUE. 
### GS:    the matrix of coefficients on the lagged control and state variables defined by output pickS and pickC.
### pickP: a matrix of coefficients that extracts a "summary of past" state vector.  Equal to ouput pickS if okPast is TRUE.
###        This "state" has all the information needed at t to predict optimally the future of the optimized system from t+1
###        onward.  May be smaller than pickS.
### okS:   input pickS works as state vector
### okC:   input pickC works with output pickS as "control" vector
### okPast: output pickS also summarizes past (but may still be redundant as a summary of the past).
### redundant: input pickS is bigger than necessary, even though okS.  Output pickS will be smaller, unless allowReduntant.
### The solution was in the form y(t)=G1 %*% y(t-1)+impact %*% z(t).  Now it's in the form pick %*% y(t)=GS %*% pick %*% y(t-1)+ PsiS %*% z(t),
### where pick=rbind(pickC,pickS).  So in the new system, control is stacked above state.  Also, new C is HS %*% new S in solution.
  REALSMALL <- 1e-9
  nr <- dim(G1)[1]
  if(dim(G1)[1]!=dim(G1)[2]) stop("G1 must be square")
  ##if nr<nc
  ##   G1=[G1;zeros(nc-nr,nc)]
  ##end
  svdGH <- svd(cbind(G1,impact))
  top <- svdGH$d > REALSMALL
  nd <- sum(top)
  us <- svdGH$u[,top,drop=FALSE]
  uu <- svdGH$u[,!top,drop=FALSE]
  vs <- svdGH$v[,top,drop=FALSE]
  vu <- svdGH$v[,!top,drop=FALSE]
  if (is.null(pickS)) {
    pickn <- t(us)
    if (dim(pickn)[1] > 0) {
      qrp <- qr(pickn)
      pickn[, qrp$pivot] <- qr.R(qrp)
    }
    ok <- FALSE
    redundant <- FALSE
  } else {
    if (is.null(dim(pickS)) ) dim(pickS) <- c(1,length(pickS))
    if (dim(pickS)[1] > dim(pickS)[2]) pickS <- t(pickS)
    svdpu <- svd(pickS %*% us)
    ok <- (sum(abs(svdpu$d) > REALSMALL) == dim(us)[2])
    ##     dp <- svdp$d[topp]
    ##     ups <- svdp$u[topp,topp,drop=FALSE]
    ##   pinw <- all( abs( pickS %*% ( diag(nr) - us %*% t(us) ) ) < REALSMALL )
    ##   winp <- all( abs( ( diag(nr) - vps %*% t(vps) ) %*% us ) < REALSMALL)  
    ##  ok <- pinw && winp
    redundant <- (dim(pickS)[1] > dim(us)[2])
    if (ok && (allowRedundant || !redundant )) {
      pickn <- pickS
    } else {
      pickn <- t(us)     
    }
  }
  pickn <- pickn / pickn[which.max(abs(pickn))] # so get positive 1's when it's possible to have all unit weights.
  if ( is.null(pickC) ) {
    picknC <- t(uu)
    if (dim(picknC)[1] > 0) {
      qruu <- qr(picknC)
      picknC[,qruu$pivot] <- qr.R(qruu)
    }
    okc <- FALSE
  } else {
    if (is.null(dim(pickC)) ) dim(pickC) <- c(1,length(pickC))
    if (dim(pickC)[1] > dim(pickC)[2]) pickC <- t(pickC)
    okc <- FALSE                     # unless conditions below are met
    if ( !is.null(pickC) && dim(pickC)[1] == nr - dim(pickn)[1] ) {
      okc <- ( sum( svd(rbind(pickC,pickn))$d > REALSMALL ) == nr ) # This checks that pickC and pickn will stack to full rank
    }
    if (okc) {
      picknC <- pickC
    } else {
      picknC <- t(uu)
      if (dim(picknC)[1] > 0) {
        qruu <- qr(picknC)
        picknC[,qruu$pivot] <- qr.R(qruu)
      }
    }
  }
  picknC <- picknC / picknC[which.max(abs(picknC))] # so get positive 1's when it's possible to have all unit weights.
  ## Check whether pickS summarizes the past
  svdG <- svd(G1)
  topG <- svdG$d > REALSMALL
  vsG <- svdG$v[, topG, drop=FALSE]
  vuG <- svdG$v[, !topG, drop=FALSE]
  pinv <- all( abs( pickn %*% vuG ) < REALSMALL )
  svdp <- svd(pickn)
  topp <- svdp$d > REALSMALL
  ndp <- sum(topp)
  vp <- svdp$v[,topp,drop=FALSE]
  vinp <- all(abs( (diag(nr) - vp %*% t(vp)) %*% vsG) < REALSMALL)
  if ( pinv && vinp ) {
    pickP <- pickn
  } else {
    pickP <- t(vsG)
    if (dim(pickP)[1] > 0) {
      qrvsG <- qr(pickP)
      pickP[, qrvsG$pivot] <- qr.R(qrvsG)
    }
  }
  pickP <- pickP / pickP[which.max(abs(pickP))]
  tmat <- rbind(picknC,pickn) 
  GS <- tmat %*% G1 %*% solve(tmat)
  PsiS <- tmat %*% impact
  HS <- t(uu) %*% solve(tmat)
  ## nh <- dim(HS)[1]
  nh <- dim(picknC)[1]
  if (nh > 0) HS <- solve(HS[1:nh, 1:nh],HS[1:nh, , drop=FALSE])
  HS <- -HS[, -(1:nh)]
  if (max(abs(HS)) > REALSMALL) HS <- HS / HS[which.max(abs(HS))]
  return(list(pickS=pickn, pickC=picknC, GS=GS, HS=HS, PsiS=PsiS, pickP=pickP,  okS=ok, okC=okc,  redundant=redundant, okPast=vinp))
}
