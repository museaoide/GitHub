setUpScrep <- function(crange=c(0.,0.5), pw, ww,
                       U=function(x,y) { log(pmax(x, .Machine$double.eps*2)) + log(pmax(y-x, .Machine$double.eps*2))}, kappa=1,
                     DxU=function(x,y)  {ifelse(x > 0 & y>x, 1/x - 1/(y-x), 0)}, lbp=1e-7, PNLTY=1) {
  ## Notes 1/3/09:  Need to distinguish systematically between "compacted" derivative terms and others.  E.g., when x <- diagmul.tensor(h, "iw", eaU),
  ## then dxdh compacted is simply eaU, whereas the full matrix version would be a 3-level tensor with zeros everywhere that iw has different values in
  ## h and eaU.  There may be places below where compacted derivatives are used as if not compacted, and probably some places where compacted derivatives
  ## could be used where they are not.
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
    ## computes objective function, including the penalty defining the constraints, and its derivative vector
    ## dc:    determines the points in crange that form the support of the c distribution.  The entries are all positive and
    ##        sum to one.  Support point i is at (sum_1^i dc_i) * (crange[2] - crange[1]) + crange[1].
    ## pc:    Marginal probabilities on the points in c's support, with last one omitted ( = 1 - sum of others).
    ## alpha: 1/(Lagrange multiplier on info constraint).  This is treated as a free parameter, so that the FOC's
    ##        let us write the joint density as pc %*% exp(alpha*U) %*% mu(w) for a pdf pc and a weighting function mu.
    require("tensorA")
    if ((sum(dc) > 1) || any( dc < 0 ) ) {
      message("bad dc vector")
      return(list(obf=-1e100))
    }
    if ((sum(pc) > 1) || any (pc < 0 ) ) {
      message("bad pc vector")
      return(list(obf=-1e100))
    }
    nc <- length(dc)
    cc <- to.tensor(cumsum(dc),dims=c(ic=nc))
    pc <- to.tensor(c(pc, 1 - sum(pc)), dims=c(ic=nc))
    cc <- cc * (crange[2] - crange[1]) + crange[1]
    ww <- to.tensor(ww,dims=c(iw=nw))
    pw <- to.tensor(pw,dims=c(iw=nw))
    ## cltw <- outer(cc, ww, function(a,b){a<b-.Machine$double.eps*2} )
    Umat <- as.tensor(outer(as.vector(cc), as.vector(ww), FUN=U), dims=c(dim(cc), dim(ww)))
    ## deriv---------------
    dUmatdcc <- as.tensor(outer(as.vector(cc),as.vector(ww), FUN=DxU), dims=c(dim(cc), dim(ww))) #compacted
    ## dUmatdcc <- to.tensor(0,c(ic=nc, iw=nw, idc=nc))
    ## for (i in 1:nc) dUmatdcc[i, , i] <- DxU(cc[i],ww)
    ## /deriv---------------
    eaU <- exp(alpha * Umat)
    ## deriv----------------
    deaUdcc <- eaU * dUmatdcc * alpha   #compacted
    ## deaUdcc <- rep(eaU, nc, pos=3, name="idc")
    ## deaUdcc <- deaUdcc * dUmatdcc * alpha
    deaUdalpha <- eaU * Umat
    ## /deriv---------------
    peaU <- diagmul.tensor(pc, "ic", eaU) 
    ## deriv----------------
    dpeaUdpc <- eaU                     #compacted
    ## dpeaUdpc <- to.tensor(0, c(ic=nc, iw=nw, idc=nc))
    ## for (i in 1:nc) dpeaUdpc[i, , i] <- eaU[i, ]
    dpeaUdcc <- diagmul.tensor(pc, "ic", deaUdcc) #compacted
    dpeaUdalpha <- peaU * Umat
    ## /deriv---------------
    pwhat <- margin.tensor(peaU, "ic")
    ## deriv----------------
    ## the first two of these, below are compacted:  cc[i] and pc[i] affect only column i of peaU
    dpwhatdcc <- dpeaUdcc               # compacted
    dpwhatdpc <- dpeaUdpc               #compacted
    dpwhatdalpha <-  margin.tensor(dpeaUdalpha,"ic")
    ## /deriv---------------
    h <- pw/pwhat                       
    ## deriv---------------
    dhdpwhat <-  -h/pwhat               #compacted.
    dhdcc <- diagmul.tensor(dhdpwhat, "iw", dpwhatdcc)
    dhdpc <- diagmul.tensor(dhdpwhat, "iw", dpwhatdpc)
    dhdalpha <- diagmul.tensor(dhdpwhat, "iw", dpwhatdalpha)
    ## /deriv--------------
    pdf <- diagmul.tensor(h, "iw", peaU)
    ## deriv--------------
    dpdfdh <- peaU  #compacted
    ## dpdfdh <- to.tensor(0, c(ic=nc, iw=nw, idw=nw))
    ## for (j in 1:nw) dpdfdh[ic=1:nc, iw=j, idw=j] <- peaU[ic=1:nc, iw=j]
    ##!! forgot the need to recognize that only diagonal elements are non-zero and messed up some lines below.
    ## dpdfdpc <- to.tensor(0,c(ic=nc,iw=nw,idc=nc))
    ## eaUh <- diagmul.tensor(h, "iw", eaU)
    dpdfdpc <- diagmul.tensor(h, "iw", dpeaUdpc)
    dpdfdpc <- diag.tensor(dpdfdpc,  by="iw", mark="d")[[icd=~idp]] + diagmul.tensor(dhdpc[[ic=~idp]], "iw", peaU) 
    ## for (i in 1:nc) dpdfdpc[ic=i,iw=1:nw , idc=1:nc] <- eaUh[ic=i, iw=1:nw]
    ## dpdfdpc <- dpdfdpc + mul.tensor(dpdfdh, "idw", dhdpc, "iw")
    dpdfdcc <- diagmul.tensor(h, "iw", dpeaUdcc)
    dpdfdcc <- diag.tensor(dpdfdcc, by="iw", mark="d")[[icd=~idc]] + diagmul.tensor(dhdcc[[ic=~idc]], "iw", peaU)
    dpdfdalpha <- diagmul.tensor(dhdalpha,"iw", peaU) + diagmul.tensor(h, "iw", dpeaUdalpha)
                                        # #t(c(dhdalpha) * t(peaU)) + t(h * t(dpeaUdalpha)) 
    ## /deriv---------------
    h0 <- margin.tensor(pdf, "ic")
    ## deriv----------------------
    dh0dcc <- margin.tensor(dpdfdcc, "ic")
    dh0dpc <- margin.tensor(dpdfdpc, "ic")
    dh0dalpha <- margin.tensor(dpdfdalpha, "ic")
    ## /deriv----------------------
    h <- pw / h0
    ## deriv----------------------
    dhdcc <- diagmul.tensor(-(h / h0 ), "iw", dh0dcc)
    dhdpc <- diagmul.tensor(-(h / h0 ), "iw", dh0dpc)
    dhdalpha <- diagmul.tensor(-(h/h0), "iw", dh0dalpha)
    ## /deriv-------------------------------
    pdf <- diagmul.tensor(h, "iw", pdf) #t(c(h) * t(pdf ))
    ## deriv--------------------------------
    dpdfdcc <- diagmul.tensor(dhdcc, "iw", pdf ) + diagmul.tensor(h, "iw", dpdfdcc)
    dpdfdpc <- diagmul.tensor(dhdpc, "iw", pdf) + diagmul.tensor(h, "iw", dpdfdpc)
    dpdfdalpha <- diagmul.tensor(dhdalpha, "iw", pdf) + diagmul.tensor(h, "iw", dpdfdalpha)
    ## /deriv--------------
    EU <- sum(pdf * Umat)
    ## deriv-----------------
    dEUdcc <- margin.tensor(diagmul.tensor(dpdfdcc, names(Umat), Umat), names(Umat))
    dEUdpc <- margin.tensor(diagmul.tensor(dpdfdpc, names(Umat), Umat), names(Umat))
    dEUdalpha <- margin.tensor(diagmul.tensor(dpdfdalpha, names(Umat), Umat), names(Umat)) #== sum(dpdfdalpha * Umat)?
    ## /deriv------------------
    ## k <- sum(pdf * log(ifelse(pdf>0,pdf,1))) - sum(pc * log(ifelse(pc>0,pc,1))) - sum(pw * log(ifelse(pw > 0, pw, 1)))
    ## paranoid version above should not be needed with lbp, and version below makes deriv easier.
    k <- sum(pdf * log(pdf)) - sum(pc * log(pc)) - sum(pw * log(pw))
    ## deriv-------------------------
    dkdcc <- margin.tensor(diagmul.tensor(dpdfdcc, names(pdf), log(pdf) ) + dpdfdcc, names(pdf))
    dkdpc <- margin.tensor(diagmul.tensor(dpdfdpc, names(pdf), log(pdf) ) + dpdfdpc, names(pdf))
    - log(pc)[[ic=~idp]] - one.tensor(c(idp=nc))
    dkdalpha <- margin.tensor(diagmul.tensor(dpdfdalpha, names(pdf), log(pdf) ) + dpdfdalpha, names(pdf))
    ## /deriv--------------------
    ## obf <- EU - lambda * k
    ## two lines below apply a penalty for pdf values near zero.  This prevents the search from wandering over large
    ## changes in large negative values of log p's, and may speed convergence.
    ##---------------------------------
    ## zerop <- pdf[ (pdf < lbp)]          #sum(pdf) is one by construction, so don't penalize its getting close to one.
    if (!any(pc < lbp)) {
      zpPnlty <- 0
      dzpPnlty <- 0
    } else {
      zerop <- pc[ pc < lbp ]
      zerop <- to.tensor(zerop, c(izp=length(zerop)))
      if (sum(pc) > 1 - lbp ) zerop <- c(zerop, 1-sum(pc))
      dzeropdp <- delta.tensor(c(idp=nc))$izp.idp
      dzeropdp <- dzeropdp[ izp= (pc < lbp), idp=1:nc, drop=FALSE]
      zpPnlty <- sum((lbp / zerop - 1)^2 ) 
      dzpPnlty <- mul.tensor(2 * (lbp / zerop -1) * (-lbp / zerop^2), "izp", dzeropdp,  "izp")
    }
    if(!any(dc < lbp)) {
      zcPnlty <- 0
      dzcPnlty <- 0
    } else {
      zerodc <- dc[dc < lbp]
      if (sum(dc) > 1 - lbp ) zerodc <- c(zerodc, 1-sum(dc))
      zerodc <- to.tensor(zerodc, c(zc=length(zerodc)))
      ## deriv------------------------
      zerodc <- to.tensor(zerodc, c(izc=length(zerodc)))
      ddcdcc <- toeplitz(nc:1)
      ddcdcc[lower.tri(ddcdcc)] <- 0
      ddcdcc <- solve(ddcdcc)
      ddcdcc <- as.tensor(ddcdcc, dims=c(iddc=nc,idc=nc))
      dzerocdcc <- delta.tensor(c(iddc=nc))$izc.iddc
      dzerocdcc <- dzerocdcc[ izc= (dc < lbp), iddc=1:nc, drop=FALSE]
      dzerocdcc <- mul.tensor(dzerocdcc, "iddc", ddcdcc, "iddc")
      dzcPnlty <- mul.tensor(2 * (lbp / zerop -1) * (-lbp / zerop^2), "izp", dzerocdcc,  "izp")
      zcPnlty <- sum((lbp / zeroc -1)^2)
    }
    ## /deriv------------------------
    ## zerodc <- dc[dc < 0]
    ## if(sum(dc) >1) zerodc <- c(zerodc, 1-sum(dc))
    ## obf <- EU  - PNLTY * (sum((lbp / zerop - 1)^2 ) + sum((lbp / zerodc -1)^2) + 100 * sum(abs(rep(1,nc) %*% pdf - pw)) + 10 * max(k-kappa,0)^2)
    ##obf <- EU  - PNLTY * (sum((lbp / zerop - 1)^2 ) + sum((lbp / zerodc -1)^2)  + 10 * max(k-kappa,0)^2)
    obf <- EU  - PNLTY * (zpPnlty + zcPnlty  + 10 * max(k-kappa,0)^2)
    ## deriv-------------------------------
    dobfdcc <- dEUdcc - PNLTY * (dzcPnlty + 10 * (k > kappa) * 2 * (k - kappa) * dkdcc)
    dobfdpc <- dEUdpc - PNLTY * (dzpPnlty + 10 * (k > kappa) * 2  * (k - kappa) * dkdpc)
    dobfdalpha <- dEUdalpha - PNLTY * (10 * (k > kappa) * 2 * (k - kappa) * dkdalpha)
    ## /deriv---------------------------------
    return(list(obf=obf, EU=EU, cc=cc, pc=pc, pdf=pdf, info=k/log(2), dobfdcc=dobfdcc, dobfdpc=dobfdpc, dobfdalpha=dobfdalpha))
  }
  screpDCO <- function(dc=rep(.1,9), pc=rep(.111,8), alpha=1) {
    ## computes objective function, including the penalty defining the constraints.
    ## dc:    determines the points in crange that form the support of the c distribution.  The entries are all positive and
    ##        sum to one.  Support point i is at (sum_1^i dc_i) * (crange[2] - crange[1]) + crange[1].
    ## pc:    Marginal probabilities on the points in c's support, with last one omitted ( = 1 - sum of others).
    ## alpha: 1/(Lagrange multiplier on info constraint).  This is treated as a free parameter, so that the FOC's
    ##        let us write the joint density as pc %*% exp(alpha*U) %*% mu(w) for a pdf pc and a weighting function mu.
    require("tensorA")
    if ((sum(dc) > 1) || any( dc < 0 ) ) {
      message("bad dc vector")
      return(list(obf=-1e100))
    }
    if ((sum(pc) > 1) || any (pc < 0 ) ) {
      message("bad pc vector")
      return(list(obf=-1e100))
    }
    nc <- length(dc)
    cc <- to.tensor(cumsum(dc),dims=c(ic=nc))
    pc <- to.tensor(c(pc, 1 - sum(pc)), dims=c(ic=nc))
    cc <- cc * (crange[2] - crange[1]) + crange[1]
    ww <- to.tensor(ww,dims=c(iw=nw))
    pw <- to.tensor(pw,dims=c(iw=nw))
    ## cltw <- outer(cc, ww, function(a,b){a<b-.Machine$double.eps*2} )
    Umat <- as.tensor(outer(as.vector(cc), as.vector(ww), FUN=U), dims=c(dim(cc), dim(ww)))
    eaU <- exp(alpha * Umat)
    peaU <-diagmul.tensor(pc, "ic", eaU) 
    pwhat <- margin.tensor(peaU, "ic")
    h <- pw/pwhat
    pdf <- diagmul.tensor(h, "iw", peaU)
    h0 <- margin.tensor(pdf, "ic")
    h <- pw / h0
    pdf <- diagmul.tensor(h, "iw", pdf) #t(c(h) * t(pdf ))
    EU <- sum(pdf * Umat)
    ## k <- sum(pdf * log(ifelse(pdf>0,pdf,1))) - sum(pc * log(ifelse(pc>0,pc,1))) - sum(pw * log(ifelse(pw > 0, pw, 1)))
    ## paranoid version above should not be needed with lbp, and version below makes deriv easier.
    k <- sum(pdf * log(pdf)) - sum(pc * log(pc)) - sum(pw * log(pw))
    ## obf <- EU - lambda * k
    ## two lines below apply a penalty for pdf values near zero.  This prevents the search from wandering over large
    ## changes in large negative values of log p's, and may speed convergence.
    ##---------------------------------
    ## zerop <- pdf[ (pdf < lbp)]          #sum(pdf) is one by construction, so don't penalize its getting close to one.
    zerop <- pc[ pc < lbp ]
    if (sum(pc) > 1 - lbp ) zerop <- c(zerop, 1-sum(pc))
    zerodc <- dc[dc < lbp]
    if (sum(dc) > 1 - lbp ) zerodc <- c(zerodc, 1-sum(dc))
    zerodc <- to.tensor(zerodc, c(zc=length(zerodc)))
    ## zerodc <- dc[dc < 0]
    ## if(sum(dc) >1) zerodc <- c(zerodc, 1-sum(dc))
    ## obf <- EU  - PNLTY * (sum((lbp / zerop - 1)^2 ) + sum((lbp / zerodc -1)^2) + 100 * sum(abs(rep(1,nc) %*% pdf - pw)) + 10 * max(k-kappa,0)^2)
    ##obf <- EU  - PNLTY * (sum((lbp / zerop - 1)^2 ) + sum((lbp / zerodc -1)^2)  + 10 * max(k-kappa,0)^2)
    obf <- EU  - PNLTY * (sum((lbp / zerop - 1)^2 ) + sum(zerodc^2)  + 10 * max(k-kappa,0)^2)
    return(list(obf=obf, EU=EU, cc=cc, pc=pc, pdf=pdf, info=k/log(2)))
  }
  return(list(dobf=screpDC, obf=screpDCO))
}
