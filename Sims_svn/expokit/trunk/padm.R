padm <- function( A, p=6 ){
  ##  PADM computes the matrix exponential exp(A) using the irreducible 
  ##  (p,p)-degree rational Pade approximation to the exponential function.
  ##
  ##  E = padm( A )
  ##
  ##  See also CHBV, EXPOKIT and the MATLAB supplied functions EXPM and EXPM1.

  ##  Roger B. Sidje (rbs@maths.uq.edu.au)
  ##  EXPOKIT: Software Package for Computing Matrix Exponentials.
  ##  ACM - Transactions On Mathematical Software, 24(1):130-156, 1998
  ##  converted to R by Chris Sims (sims@princeton.edu)

  n <- dim(A)[1]

  ## Pade coefficients (1-based instead of 0-based as in the literature)
  c <- rep(0,p+1)
  c[1]  <-  1;
  for (k in 1:p){
    c[k+1]  <-  c[k]*((p+1-k)/(k*(2*p+1-k)));
  }

  ## Scaling

  s  <-  max(apply(abs(A),1,sum));
  if (s > 0.5){ 
    s  <-  max(0,trunc(log(s)/log(2))+2)
    A  <-  2^(-s)*A
  }

  ## Horner evaluation of the irreducible fraction (see ref. above)

  I  <-  diag(1,nrow=n,ncol=n);
  A2  <-  A %*% A;
  Q  <-  c[p+1]*I;
  P  <-  c[p]*I;
  odd  <-  TRUE
  for (k in seq(p-1,1,-1)) {
    if (odd){
      Q  <-  Q %*% A2 + c[k]*I;
    } else {
      P  <-  P %*% A2 + c[k]*I;
    }
    odd  <-  !odd
  }
  if (odd){
    Q  <-  Q %*% A
    Q  <-  Q - P
    E  <-  -(I + 2*solve(Q,P))
  } else {
    P  <-  P %*% A
    Q  <-  Q - P
    E  <-  I + 2 * solve(Q,P)
  }

  ## Squaring
  if( s>1) {
    for (k in 1:s) {
      E  <-  E %*% E
    }
  }
  return(E)
}
