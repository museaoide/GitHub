sysmat <- function(B) {
## Constructs the lags*nv x lags*nv system matrix for the "stacked" first-order
## version, from the nv x nv x lags array of coefficients By returned by rfvar3.
    n <- dim(B)
    dim(B) <- c(n[1],n[2]*n[3])
    B <- rbind(B,diag(1,nrow=n[1]*(n[3]-1),ncol=n[2]*n[3]))
    return(B)
  }

