doubling <- function(A,omega,crit=1e-9) {
  V <- omega
  Aj <- A
  j <- 1
  vinc <- sum(abs(V))
  if (!is.complex(A)) {
    while (j < 10000 && vinc > crit ) {
      dv <-  Aj %*% V %*% t(Aj)
      vinc <- sum(abs(dv))
      V <- V + dv
      Aj <- Aj %*% Aj
      j <- j+1
    }
  } else {
    while (j < 10000 && vinc > crit ) {
      dv <-  Aj %*% V %*% Conj(t(Aj))
      vinc <- sum(abs(dv))
      V <- V + dv
      Aj <- Aj %*% Aj
      j <- j+1
    }
  }
  ## print(j)
  ## print(sum(abs(dv)))
  if ( vinc > crit) warning("unconverged doubling")
  return(V)
}
    
