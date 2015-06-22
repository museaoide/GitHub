kf <- function(y,H,shat,sig,G,M) {
  ## s is the state, and the plant equation is s(t)=Gs(t-1)+Me, where e is
  ## N(0,I).  The observation equation is y(t)=Hs(t).  The prior distribution for
  ## s is N(shat,sig).  To handle the standard Kalman Filter setup where the observation
  ## equation is y(t)=Hs(t)+Nu(t), expand s to become [s;v] (where v=Nu), expand H to [H I], replace
  ## G by [G 0;0 0], and replace M with [M 0;0 N].  The resulting model has no "error" in the
  ## observation equation but is equivalent to the original model with error in the state equation.
  ## The posterior on the new state is N(shatnew,signew) and lh is a two-dimensional vector containing
  ## the increments to the two component terms of the log likelihood function.  They are added 
  ## to form the log likelihood, but are used separately in constructing a concentrated or marginal
  ## likelihood. yhat is the error in the forecast of y based on the prior.lh=zeros(1,2);
  omega <- G %*% sig %*% t(G) + M %*% t(M)
  if (isTRUE(all.equal(H, 0)) ) { # No observation case.  Just propagate the state.
    lh <- c(0,0)
    shatnew <- G %*% shat
    signew <- omega
    fcsterr <- y                        # y had better be 0
    if (!all.equal(y,0) ) warning("zero H but non-zero y")
  } else {   
    svdo <- svd(omega)
    svdhud <- svd( H %*% t(t(svdo$u) * sqrt(svdo$d)) )
    if (isTRUE(all.equal(svdhud$d, 0))) { # Observation is uninformative.  Again propagate state.
      lh <- c(0,0)
      shatnew <- G %*% shat
      signew <- omega
      fcsterr <- y                      # y had better be 0
      if (!all.equal(y,0) ) warning("Uninformative H but non-zero y")
    } else {
      first0 <- match(TRUE, svdhud$d < 1e-7)
      if (is.na(first0)) first0 <- if(is.null(dim(H))) 2 else min(dim(H))+1
      u <- svdhud$u[ , 1:(first0-1), drop=FALSE]
      v <- svdhud$v[ , 1:(first0-1), drop=FALSE]
      d <- svdhud$d[1:(first0-1), drop=FALSE]
      fac <- t(t(svdo$v) * sqrt(svdo$d)) 
      fcsterr <- y-H %*% G %*% shat
      ferr <- t(t(v) * 1/d) %*% t(u) %*% fcsterr
      lh <- c(0,0)
      lh[1] <- -.5 * crossprod(ferr)
      lh[2] <- -sum( log(d) )
      shatnew <- fac %*% ferr + G %*% shat
      signew <- fac %*% (diag(dim(v)[1]) - v %*% t(v) ) %*% t(fac)
    }
  }
  return(list(shatnew=shatnew, signew=signew, lh=lh, fcsterr=fcsterr))
}
