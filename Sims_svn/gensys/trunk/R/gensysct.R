gensysct <- function(g0, g1, c0=rep(0,dim(g0)[1]), psi, pi, div=-1) {
  ##System given as
  ##        g0*Dy(t)=g1*y(t)+c+psi*epsilon(t)+pi*eta(t),
  ##with epsilon an exogenous variable process and eta being endogenously determined
  ##white noise expectational errors.  Returned system is
  ##       Dy(t) =  -G1  %*% y(t) +impact %*% epsilon(t) + C ,
  ## for the case where epsilon is white noise.  If epsilon is not white noise, the solution
  ## is more complicated but can be computed from the returned qz decomposition.
  ## If div is omitted from argument list, a div>0 is calculated.
  ## Also returned is the qz decomposition, q'az'=g0, q'bz'=g1, with a and b
  ## upper triangular and the system ordered so that all zeros on the diagonal of b are in
  ## the lower right corner, all cases where the real part of bii/aii is greater than or 
  ## equal to div appear in the next block above the zeros, and the remaining bii/aii's 
  ## all have bii/aii<div .    See the paper "Solving Linear Rational Expectations Models", 
  ## http://eco-072399b.princeton.edu/yftp/gensys .  Note that if one simply wants the backward
  ## and forward projection of y on epsilon, ignoring existence and uniqueness questions, the
  ## projection can be computed by Fourier methods.
  ## eu[1]=(existence); eu[2]=(uniqueness). Else eu==c(-2,-2)  => indeterminacy via singularity
  ## in the equation system.
  ## derivs:  If >=0, the forward part of the solution involves delta functions (derivs=0) or derivatives
  ## up to order derivs, of the expected future path of the exogenous process epsilon.  The paper cited above
  ## explains how to derive the forward solution from this function's return values.
  vnames <- dimnames(g0)[[2]]
  shocknames <- dimnames(psi)[[2]]
  experrnames <- dimnames(pi)[[2]]
  realsmall <- 1e-7
  eu <- c(0,0)
  fixdiv <- (div > 0)
  n <- dim(g0)[1]
  if ( !is.null(dim(psi)) ) {
    nshock <- dim(psi)[2]
  } else {
    if (is.null(psi)) {
      nshock <- 0
    } else {
      nshock <- 1
    }
  }
  qzout <- qz(g0,g1)
  zxz <- any((abs(diag(qzout$a)) < realsmall) & (abs(diag(qzout$b)) < realsmall))
  if (zxz) {
    "Coincident zeros.  Indeterminacy and/or nonexistence.\n"
    eu <- c(-2,-2)
    return(list(eu=eu,qzdec=qzout))
  }
  zeroax <- abs(diag(qzout$a)) < realsmall
  nzero <- sum(zeroax)
  unstabx <- cos(Arg(diag(qzout$b))-Arg(diag(qzout$a))) > realsmall # near zero roots don't count
  unstabx <- (!zeroax) & unstabx
  if (!fixdiv) {
    if (! any(unstabx)) {
      div <- .001
    } else {
      div <- .5 * min(Re(diag(qzout$b)[unstabx] / diag(qzout$a)[unstabx]))
    }
  }
  ## Now that we know what div is, and are sure of no double zeros, reset unstabx
  unstabx <- zeroax | (Re(diag(qzout$b) / (zeroax+diag(qzout$a))) > div)
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
  q1 <- qq[six, , drop=FALSE]
  q2 <- qq[uix, , drop=FALSE]
  z1 <- t(Conj(qzout$z[, six, drop=FALSE]))
  z2 <- t(Conj(qzout$z[, uix, drop=FALSE]))
  a2 <- qzout$a[uix,uix,drop=FALSE]
  b2 <- qzout$b[uix,uix,drop=FALSE]
  etawt <- q2 %*% pi
  neta <- if (is.matrix(pi)) dim(pi)[2] else if (is.null(pi)) 0 else 1
  ndeta <- min(nunstab,neta)
  if(ndeta==0){
    ueta <- matrix(0,nunstab,0)
    deta <- vector("numeric",0)
    veta <- matrix(0,neta,0)
    bigev <- vector("logical",0)
  } else {
    sd <- svd(etawt)
    ueta <- sd$u; deta <- sd$d; veta <- sd$v
    bigev <- deta>realsmall
    ueta<-ueta[,bigev,drop=FALSE]
    veta<-veta[,bigev,drop=FALSE]
    deta<-deta[bigev]
  }
  eu[1] <- sum(bigev) >= nunstab
  ##----------------------------------------------------
  ## Note that existence and uniqueness are not just matters of comparing
  ## numbers of roots and numbers of endogenous errors.  These counts are
  ## reported below because usually they point to the source of the problem.
  ##------------------------------------------------------
  etawt1 <- q1 %*% pi
  ndeta1 <- min(n-nunstab,neta)
  if(ndeta1==0){
    ueta1 <- matrix(0,n-nunstab,0)
    deta1 <- vector("numeric",0)
    veta1 <- matrix(0,neta,0)
    bigev1 <- vector("logical",0)
  } else {
    sd <- svd(etawt1)
    ueta1<-sd$u
    deta1 <- sd$d
    veta1 <- sd$v
    bigev1 <- deta1 > realsmall
  }
  if (any(bigev1)) { #needed because empty dimensions are dropped after select
    ueta1 <- ueta1[,bigev1,drop=FALSE]
    veta1 <- veta1[,bigev1,drop=FALSE]
    deta1 <- deta1[bigev1]
    loose <- veta1-veta %*% t(Conj(veta)) %*% veta1
    svdl <- svd(loose)
    loose <- sum(abs(svdl$d)>realsmall*n)
    unq <- (loose==0)
  } else {
    ueta1 <- matrix(1,n-nunstab,0)
    veta1 <- matrix(1,neta,0)
    deta1 <- vector("complex",0)
    unq <- TRUE
  }
  if (unq) {
    eu[2] <- 1
  } else
  {
    cat("Indeterminacy.", loose, "loose endog errors.\n")
  }
  ## Note: if v is a vector of length n and m is an nxp matrix,
  ## v*m==diag(v)%*%m, m/v==solve(diag(v),m)==diag(v)\m (matlab notation)
  ##
  tmat <- cbind(diag(n-nunstab),
                -t(Conj((ueta %*% (t(Conj(veta))/deta)) %*% veta1 %*% (deta1 * t(Conj(ueta1)))))  )  
  G0<- rbind( tmat %*% qzout$a, cbind(matrix(0,nunstab,n-nunstab), diag(nunstab)))
  G1<- rbind(tmat %*% qzout$b, matrix(0,nunstab,n))
  ##----------------------
  ## G0 is always non-singular because by construction there are no zeros on
  ## the diagonal of a[1:(n-nunstab),1:(n-nunstab)], which forms G0's ul corner.
  ##-----------------------
  G0I <- solve(G0)
  G1 <- G0I%*%G1
  ##----------- uix can be empty, e.g. in indeterminate systems with no unstable roots ------------
  if(is.null(uix)){
    C <- G0I %*% tmat %*% qq %*% c0
    ## fmat <- matrix(0,0,0)
    ## fwt <- matrix(0, 0, nshock)
    impact <- G0I %*% tmat %*% qq %*% psi
  }else {
    ## line below is for the discrete time case
    ## C <- G0I %*% rbind(tmat%*% qq %*%c0,solve(qzout$a[uix,uix,drop=FALSE]-qzout$b[uix,uix,drop=FALSE],q2%*%c0) )
#    C <- G0I %*% rbind(tmat%*% qq %*%c0,solve(-qzout$b[uix,uix,drop=FALSE],q2%*%c0) )
    C <- G0I %*% rbind(tmat%*% qq %*%c0, matrix(0,nunstab,1) ) # WYP, 8/29/2008, plug in zeros for the unstable root part
    impact <- G0I %*% rbind(tmat %*% qq %*% psi, matrix(0,nunstab, nshock))
    ## fmat <- solve(qzout$b[uix,uix,drop=FALSE],qzout$a[uix,uix,drop=FALSE])
    ## fwt <- -solve(qzout$b[uix,uix,drop=FALSE],q2 %*% psi)
  }
  ## ywt <- G0I[,uix,drop=FALSE]
  loose <- etawt1 %*% (diag(neta) - veta %*% t(Conj(veta))) 
  ##-------------------- above are output for system in terms of z'y -------
  G1<-Re(qzout$z %*% G1 %*% t(Conj(qzout$z)))
  C <- Re(qzout$z%*%C)
  impact <- Re(qzout$z%*%impact)
  ## ywt <- qzout$z%*%ywt
  loose <- Re(qzout$z %*% rbind(loose,matrix(0,nunstab,neta)))
  ## check for whether derivatives are present in forward solution
  ordera <- -1
  if (nzero > 0) {
    a33 <- qzout$a[zeroax, zeroax, drop=FALSE]
    diag(a33) <- 0 # This should already be almost exactly zero.  Prevents an infinite loop below from rounding error.
    ordera <- 0
    an <- a33
    while (sum(abs(an)) > sqrt(realsmall * nzero^2)) {
      an <- a33 %*% an
      ordera <- ordera+1
    }
  }
  dimnames(G1) <- list(vnames,vnames)
  dimnames(C) <- list(vnames, "const")
  dimnames(impact) <- list(vnames,shocknames)
  dimnames(loose) <- list(vnames,NULL)
  return(list(G1=G1, C=C, impact=impact, qzdec=qzout, eu=eu, loose=loose, derivs=ordera))
}
