Mtrop <- function(lhfcn, x0, drawjump, nit, ...) {
  mcout <- matrix(0, nit, length(x0))
  mcout[1, ] <- x0
  lhlist <- vector("numeric", nit)
  lh0 <- lhfcn(x0, ...)
  if (attr(lh0, "bad")) return("bad initial parameter")
  lhlist[1] <- lh0
  for ( it in 2:nit) {
    x1 <- drawjump(x0)
    lh1 <- lhfcn(x1, ...)               #remember lhfcn returns *minus* log lh
    if ((lh0 - lh1 > log(runif(1))) && !attr(lh1, "bad")) {
      lhlist[it] <- lh1
      lh0 <- lh1
      x0 <- x1
      mcout[it, ] <- x1
    } else {
      lhlist[it] <-  lh0
      mcout[it, ] <- x0
    }
  }
  return(list(mcout=mcout, lhlist=lhlist))
}
