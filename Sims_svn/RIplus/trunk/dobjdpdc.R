dobjdpdc <- function(p, cc, g, ww, U=function(x,y) {log(ifelse(x>0,x,1e-100)) + log(ifelse(y > x, y - x, 1e-100))}, lambda=1,
                     cderiv=TRUE, DxU=function(x,y) {1/x - 1/(y-x)}) {
  Umat <- outer(cc,ww,U)
  eU <- exp(Umat)
  peU <- p %*% eU
  h <- g/peU
  h <- c(h)
  dhdp <- - t(eU %*% diag(c(h/peU)))
  eUU <- eU * Umat
  dobjdpv <- eUU %*% h + t(p %*% eUU %*% dhdp)
  dobjdpv <-  dobjdpv - lambda * (diag(log(p)) %*% eU %*% h + eUU %*% h + eU %*% (h * log(h)) + eU %*% h
                                  + t((p * log(p)) %*% eU %*% dhdp + p %*% eUU %*% dhdp + p %*% eU %*% diag(log(h)) %*% dhdp)
                                  + t(p %*% eU %*% dhdp) -1 - log(p))
  if (cderiv) {                       #Get deriv w.r.t. cc points also
    DcU <- outer(cc,ww,DxU)
    pdeuh <- diag(p) %*% (DcU * eU) %*% diag(h)
    peuh <- diag(p) %*% eU %*% diag(h)
    dobjdcv <- diag(p) %*% DcU %*% h - lambda * ( apply(log(peuh) * pdeuh, 1, sum)  - pdeuh %*% log(apply(peuh,2,sum)) )
    return(list(dp=dobjdpv, dc=dobjdcv))
  } else {
    return(dobjdpv)
  }
}
