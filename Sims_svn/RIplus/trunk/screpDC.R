setUpScrep <- function(crange=c(0.,0.5), nc,  pw, ww, U=function(x,y) { log(pmax(x, .Machine$double.eps*2)) + log(pmax(y-x, .Machine$double.eps*2))}, kappa=1,
                     DxU=function(x,y)  {ifelse(x > 0 & y>x, 1/x - 1/(y-x), 0)}, lbp=1e-7, PNLTY=1) {
  ## Notes 5/28/2007:  derivative is not working?  needs to be made consistent in obf part with level fcn.  pdf emerges from csminwel not respecting c<w in first column.
  ## Notes 6/1/2007:  Still need to fix derivative.  With numgrad seems possibly to be working, but numerically rough.
  ## ----------------FIX ABOVE---------------------
  ## ww,pw: Points of support and marginal probabilities for w distribution.  
  ## U:     objective function U(x,w). 
  ## DxU:   derivative of U w.r.t. its first argument. If this is NULL, numerical gradient is used.
  ## crange: range of values for x 
  ## nc:     number of distinct c values.  (length of dc)
  ## lambda: multiplier on the information constraint, utility cost of capacity
  ## lbp:    lower bound on elements of dc, pc, dw, pw and 1 - upper bound on their sums.
  ## PNLTY:  weight on penalty function keeping solution away from boundary.
  if (any(diff(ww) <= 0) ) stop("ww must be increasing")
  if (any(pw < 0) || sum(pw) >1) stop("pw must be a probability vector")
  nw <- length(ww)
  ## cderiv <- !is.null(DxU)
  screpDC <- function(dc=rep(.1,9), pc=rep(.111,8), alpha=1) {
    ## computes objective function E[U] minus lambda * (information constraint) - penalty for closeness to boundaries
    ## dc:    determines the points in crange that form the support of the c distribution.  The entries are all positive and
    ##        sum to one.  Support point i is at (sum_1^i dc_i) * crange.
    ## pc:    Marginal probabilities on the points in c's support, with last one omitted (=1 - sum of others).
    ## alpha: 1/(Lagrange multiplier on info constraint).  This is treated as a free parameter, so that the FOC's
    ##        let us write the joint density as pc %*% exp(alpha*U) %*% mu(w) for a pdf pc and a weighting function mu.
    if ((sum(dc) > 1) || any( dc < 0 ) ) {
      message("bad dc vector")
      return(list(obf=-1e100))
    }
    if ((sum(pc) > 1) || any (pc < 0 ) ) {
      message("bad pc vector")
      return(list(obf=-1e100))
    }
    cc <- cumsum(dc)
    pc <- c(pc, 1 - sum(pc))
    cc <- cc * (crange[2] - crange[1]) + crange[1]
    cltw <- outer(cc, ww, function(a,b){a<b-.Machine$double.eps*2} )
    Umat <- outer(cc,ww, FUN=U)
    pdf <- exp(Umat * alpha)
    pdf <- pc * pdf
    h <- rep(1,nc) %*% pdf
    h <- pw / h
    pdf <- t(c(h) * t(pdf ))
    EU <- sum(pdf * Umat)
    k <- sum(pdf * log(ifelse(pdf>0,pdf,1))) - sum(pc * log(ifelse(pc>0,pc,1))) - sum(pw * log(ifelse(pw > 0, pw, 1)))
    ## obf <- EU - lambda * k
    ## two lines below apply a penalty for pdf values near zero.  This prevents the search from wandering over large
    ## changes in large negative values of log p's, and may speed convergence.
    ##---------------------------------
    zerop <- pdf[ cltw & (pdf < lbp)]
    ## zerodc <- dc[dc < lbp]
    ## if (sum(dc) > 1 - lbp ) zerodc <- c(zerodc, 1-sum(dc))
    zerodc <- dc[dc<0]
    if(sum(dc) >1) zerodc <- c(zerodc, 1-sum(dc))
    ## obf <- EU  - PNLTY * (sum((lbp / zerop - 1)^2 ) + sum((lbp / zerodc -1)^2) + 100 * sum(abs(rep(1,nc) %*% pdf - pw)) + 10 * max(k-kappa,0)^2)
    ##obf <- EU  - PNLTY * (sum((lbp / zerop - 1)^2 ) + sum((lbp / zerodc -1)^2)  + 10 * max(k-kappa,0)^2)
    obf <- EU  - PNLTY * (sum((lbp / zerop - 1)^2 ) + sum(zerodc^2)  + 10 * max(k-kappa,0)^2)
    ##---------------------------------
    return(list(obf=obf, EU=EU, cc=cc, pc=pc, pdf=pdf, info=k/log(2)))
  }

  dobjdpdc <- function(dc=rep(.1,9), pc=rep(.111,8)) {
    ## returns derivative and function value
    ## computes objective function E[U] minus lambda * (information constraint) - penalty for closeness to boundaries
    ## dc:    determines the points in crange that form the support of the c distribution.  The entries are all positive and
    ##        sum to one.  Support point i is at (sum_1^i dc_i) * crange.
    ## pc:    Marginal probabilities on the points in c's support, with last one omitted (=1 - sum of others).
    ##
    if ((sum(dc) > 1) || any( dc < 0 ) ) {
      message("bad dc vector")
      return(list(obf=-1e100, grad=rep(0,length(dc+length(pc)))))
    }
    if ((sum(pc) > 1) || any (pc < 0 ) ) {
      message("bad pc vector")
      return(list(obf=-1e100, grad=rep(0,length(dc+length(pc)))))
    }
    cc <- cumsum(dc)
    pc <- c(pc, 1 - sum(pc))
    cc <- cc * (crange[2] - crange[1]) + crange[1]
    cltw <- outer(cc, ww, function(a,b){a<b-.Machine$double.eps*2} )
    Umat <- outer(cc,ww, FUN=U)
    eU <- exp(Umat)
    eUU <- eU * Umat
    peU <- pc * eU
    h <- rep(1,nc) %*% peU
    h <- pw / h
    dhdp <- -t(eU %*% diag(c(h/peU))) #so w indexes rows, c indexes columns on this 
    pdf <- t(c(h) * t(peU))
    EU <- sum(pdf * Umat)
    k <- sum(pdf * log(ifelse(pdf>0,pdf,1))) - sum(pc * log(ifelse(pc>0,pc,1))) - sum(pw * log(ifelse(pw > 0, pw, 1)))
    ## obf <- EU - lambda * k
    ## two lines below apply a penalty for pdf values near zero.  This prevents the search from wandering over large
    ## changes in large negative values of log p's, and may speed convergence.
    ##---------------------------------
    zerop <- pdf[ cltw & (pdf < lbp)]
    ## zerodc <- dc[dc < lbp]
    ## if (sum(dc) > 1 - lbp ) zerodc <- c(zerodc, 1-sum(dc))
    zerodc <- dc[dc<0]
    if(sum(dc) >1) zerodc <- c(zerodc, 1-sum(dc))
    ## obf <- EU  - PNLTY * (sum((lbp / zerop - 1)^2 ) + sum((lbp / zerodc -1)^2) + 100 * sum(abs(rep(1,nc) %*% pdf - pw)) + 10 * max(k-kappa,0)^2)
    ##obf <- EU  - PNLTY * (sum((lbp / zerop - 1)^2 ) + sum((lbp / zerodc -1)^2)  + 10 * max(k-kappa,0)^2)
    obf <- EU  - PNLTY * (sum((lbp / zerop - 1)^2 ) + sum(zerodc^2)  + 10 * max(k-kappa,0)^2)
    ##---------------------------------    cc <- cumsum(dc)
    dobjdpv <- eUU %*% h + t(pc %*% eUU %*% dhdp)
    dobjdpv <-  dobjdpv - PNLTY * (k > kappa) * (diag(log(pc)) %*% eU %*% h + eUU %*% h + eU %*% (h * log(h)) + eU %*% h
                                    + t((pc * log(pc)) %*% eU %*% dhdp + pc %*% eUU %*% dhdp + pc %*% eU %*% diag(log(h)) %*% dhdp)
                                    + t(pc %*% eU %*% dhdp) -1 - log(pc))
    dobjdpv <- dobjdpv - PNLTY *  ( (pc < lbp) * (-2) * (lbp/p - 1) * lbp / p^2  +  2*(dc<0)*dc )
    DcU <- outer(cc,ww,DxU)
    pdeuh <- diag(pc) %*% (DcU * eU) %*% diag(h)
    peuh <- diag(pc) %*% eU %*% diag(h)
    dobjdcv <- diag(pc) %*% DcU %*% h - PNLTY * (k > kappa) * ( apply(plogp(peuh), 1, sum)  - pdeuh %*% log(apply(peuh,2,sum)) )
    dobjdcv <- dobjdcv - PNLTY * (dc < lbp) * (-4) * (lbp/dc - 1)^3 * lbp / dc^2
    return(list(dp=dobjdpv, dc=dobjdcv, obf=obf, EU=EU, cc=cc, pc=pc, pdf=pdf, info=k/log(2)))
  }
  screpDCS <- function(x) {
    dc <- x[1:nc]
    pc <- x[-(1:nc)]
    return(-screpDC(dc,pc)$obf)
  }
  DscrepDCS <- function(x) {
    dc <- x[1:nc]
    pc <- x[-(1:nc)]
    return(list(g=c(-dobjdpdc(dc,pc)$dc, -dobjdpdc(dc,pc)$dp[-nc]), badg=FALSE))
  }
  ## since dobjdpdc also computes the function value, it would be more efficient to have csminwel
  ## recognize this.  But since csminwel evaluates the gradient once, and the function often many times, on each
  ## iteration, the efficiency gain might not be huge.
  return(list(screpDC=screpDC, DscrepDC=dobjdpdc, screpDCS=screpDCS, DscrepDCS=DscrepDCS))
}
plogp <- function(p) {p * log( ifelse(p>0, p,1)) }
