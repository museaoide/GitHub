Pnlty <- function(g0, g1, c0=rep(0,dim(g0)[1]), psi, pi, div){
# returns a penalty for explosive roots which is increasing in the
# size of the explosive roots made stable by pushing up div (0).
# All from gensysct.R

realsmall <- 1e-7
#A <- 1e+10
A <- 1e+12

  qzout <- qz(g0,g1)
  zeroax <- abs(diag(qzout$a)) < realsmall
  nzero <- sum(zeroax)
  unstabx <- cos(Arg(diag(qzout$b))-Arg(diag(qzout$a))) > realsmall # near zero roots don't count
  unstabx <- (!zeroax) & unstabx

  n <- dim(g0)[1]
  fixdiv <- (div > 0)
  if (!fixdiv) {
    if (! any(unstabx)) {
      div <- .001
    } else {
      div <- .5 * min(Re(diag(qzout$b)[unstabx] / diag(qzout$a)[unstabx]))
    }
  }
    
roots <- Re(diag(qzout$b) / (zeroax+diag(qzout$a)))

  ## Now that we know what div is, and are sure of no double zeros, reset unstabx
  unstabx <- zeroax | (roots > div)
  nunstab <- sum(unstabx)
  ## Note that qzdivct first puts all singularities in a in lower right, then puts unstable
  ## roots on top of those.
  qzout <- qzdivct(div,qzout)
  qq <- t(Conj(qzout$q))                #to match matlab convention
  if (nunstab==n){
    six <- NULL
    uix <- 1:n
  } else
  {
    if (nunstab==0){
      uix <- NULL
      six <- 1:n
    } else
    {
      uix <- (n-nunstab+1):n
      six <- 1:(n-nunstab)
    }
  }
  #q1 <- qq[six, , drop=FALSE]
  q2 <- qq[uix, , drop=FALSE]
  #z1 <- t(Conj(qzout$z[, six, drop=FALSE]))
  #z2 <- t(Conj(qzout$z[, uix, drop=FALSE]))
  #a2 <- qzout$a[uix,uix,drop=FALSE]
  #b2 <- qzout$b[uix,uix,drop=FALSE]
  etawt <- q2 %*% pi
  neta <- if (is.matrix(pi)) dim(pi)[2] else if (is.null(pi)) 0 else 1
  ndeta <- min(nunstab,neta)
  if(ndeta==0){
    #ueta <- matrix(0,nunstab,0)
    #deta <- vector("numeric",0)
    #veta <- matrix(0,neta,0)
    bigev <- vector("logical",0)
  } else {
    sd <- svd(etawt)
    ueta <- sd$u; deta <- sd$d; veta <- sd$v
    bigev <- deta>realsmall
    #ueta<-ueta[,bigev,drop=FALSE]
    #veta<-veta[,bigev,drop=FALSE]
    #deta<-deta[bigev]
  }
nbigev <- sum(bigev)
newdiv <- div
pnlty <- 0
proots <- sort(roots[roots>div])
if ((nunstab==nbigev+2)&&(proots[1]-proots[2]+1==1)) {
	newdiv <- proots[1] + .Machine$double.eps*100
	pnlty <- A*2*(proots[1]-div)^2
} else {
	if (nunstab==nbigev+1) {
		newdiv <- proots[1] + .Machine$double.eps*100
		pnlty <- A*(proots[1]-div)^2
	}
}

return(list(newdiv=newdiv,pnlty=pnlty,nunstab=nunstab))
}