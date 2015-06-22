penalty <- function(value, boundary, threshold, power=2, cap=1e20) {
  ## function is zero, with zero derivative, at threshhold or value on same side of
  ## boundary as threshold but farther away.  It goes to cap continuously as
  ## value goes to the boundary.  As size increases it gets flatter at threshhold,
  ## rises more steeply at boundary (becomes more L-shaped).  For use with
  ## routines that use second-order information, should take power >= 2.
  bstar <- boundary + .01*(boundary - threshold)
  y <- (threshold - value)/(value - bstar)
  yb <- (threshold - boundary)/(boundary - bstar)
  pval <- ifelse(sign(value - threshold) == sign(boundary - threshold), 0,  min((y/yb)^power * cap, cap))
  ## something still fishy
  return(pval)
}
