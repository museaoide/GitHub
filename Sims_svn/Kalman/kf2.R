kf2 <- function(y,H,shat,sig,G,M) {
  ## s is the state, and the plant equation is s(t)=G %*% s(t-1)+t(M) %*% e, where e is
  ## N(0,I).  The observation equation is y(t)=Hs(t).  The prior distribution for
  ## s is N(shat,sig).  To handle the standard Kalman Filter setup where the observation
  ## equation is y(t)=Hs(t)+Nu(t), expand s to become [s;v] (where v=Nu), expand H to [H I], replace
  ## G by [G 0;0 0], and replace M with [M 0;0 N].  The resulting model has no "error" in the
  ## observation equation but is equivalent to the original model with error in the state equation.
  ## The posterior on the new state is N(shatnew,signew) and lh is a two-dimensional vector containing
  ## the increments to the two component terms of the log likelihood function.  They are added 
  ## to form the log likelihood, but are used separately in constructing a concentrated or marginal
  ## likelihood. fcsterr is the error in the forecast of y based on the prior.
  SMALLSV <- 1e-7
  omega <- G %*% sig %*% t(G) + crossprod(M)
  if (is.null(dim(H))) dim(H) <- c(1,length(H))
  nobs <- dim(H)[1]
  nstate <- dim(H)[2]
  ## stopifnot (nstate >= nobs)
  ##------------ Don't need separate treatment of H == 0.  H %*% G %*% t(H) = 0 covers it.
  ##   if (isTRUE(all.equal(H, 0))) { # No observation case.  Just propagate the state.
  ##     lh <- c(0,0)
  ##     shatnew <- G %*% shat
  ##     signew <- omega
  ##     fcsterr <- y                        # y had better be 0
  ##     if (!all.equal(y,0) ) warning("zero H but non-zero y")
  ##   } else {
  ho <- H %*% omega
  svdhoh <- svd( ho %*% t(H) )
  if (all(svdhoh$d < SMALLSV)) { # Observation is uninformative. Propagate state.
    lh <- c(0,0)
    shatnew <- G %*% shat
    signew <- omega
    fcsterr <- y - H %*% G %*% shat     # had better be 0
    if (!all(abs(fcsterr) < 1e-7)) warning("Uninformative H but non-zero fcsterr")
  } else {
    first0 <- match(TRUE, svdhoh$d < SMALLSV)
    if (is.na(first0)) first0 <- min(dim(H))+1
    u <- svdhoh$u[ , 1:(first0-1), drop=FALSE]
    v <- svdhoh$v[ , 1:(first0-1), drop=FALSE]
    d <- svdhoh$d[1:(first0-1), drop=FALSE]
    fcsterr <- y-H %*% G %*% shat
    hohifac <- (1/sqrt(d)) * t(u)       #diag(1/sqrt(d)) %*% t(u)
    ferr <- hohifac %*% fcsterr
    lh <- c(0,0)
    lh[1] <- -.5 * crossprod(ferr)
    lh[2] <- -.5 * sum( log(d) ) - .5 * log(2 * pi)
    ## log(2 * pi) term added 2013.12.11.  
    hohoifac <-hohifac %*% ho 
    shatnew <- crossprod(hohoifac, ferr) + G %*% shat
    signew <- omega - crossprod(hohoifac)
  }
  return(list(shat=shatnew, sig=signew, lh=lh, fcsterr=fcsterr))
}
