screpDCS <- function(x, ...){
  ## a wrapper for screpDC(), to let it work as the function argument to csminwel.
  nc <- (length(x)+1)/2
  dc <- x[1:nc]
  pc <- x[nc + (1:(nc-1))]
  value=-screpDC(dc=dc, pc=pc, ...)$obf; return(ifelse(is.finite(value),value,1e20))
}

DscrepDCS <- function(x, ...) {
  nc <- (length(x)+1)/2
  cc <- cumsum(x[1:nc])
  pc <- x[nc+(1:(nc-1))]
  pc <- c(pc,1-sum(pc))
  dout <- dobjdpc(p=pc, cc=cc )
  dp <- dout$dp[1:(nc-1)]-dout$dp[nc]
  
