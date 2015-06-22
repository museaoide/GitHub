rmpyA0 <- function(parvec) {
## identification scheme is
  ## xx00
  ## xxxx
  ## 00xx
  ## 000x
  ## equations are rows:  Fed policy, MD/financial, 2x2 normalized "sluggish" py sector.
  ## pindex <- c(1,2,5,6,10,11,14:16) # full set of overid'ing restrictions
  pindex <- c(1,2,5,6,10,11,13:16) # allowing contemporaneous ip in Fed equation
  A0 <- rep(0,16)
  A0[pindex] <- parvec
  A0 <- matrix(A0,4,4)
  return(A0)
}
