ksmooth <- function(btt,Stt,bnT,SnT,A,omega) {
  ## Smoothing recursion.  State evolution equation is
  ##    bn=A*bt+e,  Var(e)=omega
  ##    bt|t ~ N(btt,Stt) -- from Kalman Filter
  ##    bn|T ~ N(bnT,SnT) -- distribution of bn given full sample. From 
  ##                         KF if n=T, otherwise from this recursion
  ##    bt|T ~ N(btT,StT)
  AS=A %*% Stt;
  G <- AS %*% t(A) + omega
  ##****************************
  ##SAGI <- AS'/G
  ## line below may slightly slow the routine, but makes it robust vs.
  ## cases where part or all of the state vector is known with certainty.
  svdg <- svd(G)
  di <- svdg$d
  ##nzndx <- di > sqrt(.Machine$double.eps)
  nzndx <- di > .Machine$double.eps
  di[nzndx] <- 1/di[nzndx]
  di[-nzndx] <- 0
  SAGI <- t(AS) %*% svdg$u %*% diag(di) %*% t(svdg$v)
  ##SAGI <- t(qr.solve(t(G),AS))
  ##*****************************
  btT <- SAGI %*% (bnT - A %*% btt) + btt
  StT <- Stt - SAGI %*% AS + SAGI %*% SnT %*% t(SAGI)
  return(list(btT=btT, StT=StT))
}
