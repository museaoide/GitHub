schur <- function(A) {
  ## Returns Q and T, T upper triangular, t(Conj(Q)) %*% Q = I, Q %*% T %*% t(Conj(Q)) = A.
  ## code could be more transparent and maybe more efficient if it called zgees.f, the lapack general Schur routine,
  ## directly.  
  if (is.null(dim(A))) dim(A) <- c(1,1)
  if ( dim(A)[1] ==0) return(list(X=matrix(0,0,0), info=0))
  if ( !is.loaded("zhseqr")) dyn.load("/usr/lib/liblapack.so")
  n <- dim(A)[1]
  hessA <- .Fortran("zgehrd", as.integer(n), as.integer(1), as.integer(n),
                    A=as.complex(A), as.integer(n), tau = complex(n - 1), work = complex( n * n ),
                    lwork = as.integer(n*n), info=integer(1))
  QA <- .Fortran("zunghr", as.integer(n), as.integer(1), as.integer(n), Q = as.complex(hessA$A), as.integer(n),
                 as.complex(hessA$tau), work=complex(n * n), lwork = as.integer(n * n), info=integer(1))
  schurA <- .Fortran("zhseqr", "S", "V", as.integer(n), as.integer(1), as.integer(n), T=as.complex(hessA$A),
                     as.integer(n), w=complex(n), q=as.complex(QA$Q), as.integer(n),
                     work=complex(6*n), lwork=as.integer(6*n), info=integer(1))
  Q <- matrix(schurA$q, n)
  T <- matrix(schurA$T, n)
  return(list(Q=Q, T=T, info=schurA$info))
}
         
