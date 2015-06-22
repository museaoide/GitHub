cdfwtd <- function(weights, vals) {
  so <- order(vals)
  cdf <- cbind(vals[so], cumsum(weights[so]))
}