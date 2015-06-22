Ucapm <- function(a,b) {
  gamma <- .001
  xy <- crossprod(a,b)
  u <- xy - .5 * xy^2 - .5 * gamma * crossprod(a[-1])
  ## gamma is interpretable as the variance of an irreducible, zero-mean source of uncertainty
  ## in all assets except the first (risk free) asset.
  ## first version of this function applied it to all assets, and forgot the .5.
}