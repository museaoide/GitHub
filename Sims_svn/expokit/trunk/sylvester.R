sylvester <- function(A, C) {
  ## solves A %*% X + X %*% t(A) = C  for X.  A and C must be square.
  if(is.null(dim(A))) dim(A) <- c(1,1)
  if ( dim(A)[1] == 0) return(list(X=matrix(0,0,0), info=0))
  if(!is.loaded("ztrsyl")){
    dyn.load("/usr/lib/liblapack.so", now=FALSE)
  }
  n <- dim(A)[1]
  hessA <- .Fortran("zgehrd", as.integer(n), as.integer(1), as.integer(n),
                    A=as.complex(A), as.integer(n), tau = complex(n - 1), work = complex( n * n ),
                    lwork = as.integer(n*n), info=integer(1))
  QA <- .Fortran("zunghr", as.integer(n), as.integer(1), as.integer(n), Q = as.complex(hessA$A), as.integer(n),
                 as.complex(hessA$tau), work=complex(n * n), lwork = as.integer(n * n), info=integer(1))
  schurA <- .Fortran("zhseqr", "S", "V", as.integer(n), as.integer(1), as.integer(n), T=as.complex(hessA$A),
                     as.integer(n), w=complex(n), q=as.complex(QA$Q), as.integer(n),
                     work=complex(6*n), lwork=as.integer(6*n), info=integer(1))
  q <- matrix(schurA$q, n)
  Ctran <- t(Conj(q)) %*% C %*% q
  sylo <- .Fortran("ztrsyl", "N", "C", as.integer(1), as.integer(n), as.integer(n), as.complex(schurA$T),
                   as.integer(n), as.complex(schurA$T), as.integer(n), X=as.complex(Ctran), as.integer(n),
                   scale=double(1), info=integer(1))
  X <- matrix(sylo$X, n)/sylo$scale
  return( list( X = q %*% X %*% t(Conj(q)), info = sylo$info ))
}
