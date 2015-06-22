Umix <- function(a, b) {
  vm <- a - b - 1
  vp <- a - b + 1
  dnm <- dnorm(vm)
  dnp <- dnorm(vp)
  uve <- dnm + dnp
  uv <- log(uve)
  ddx <- (1/uve) * (-vm * dnm - vp * dnp)
  attr(uv, "gradient") <- ddx
  return(uv)
}