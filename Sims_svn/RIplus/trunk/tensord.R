tensord <- function(A,B,alongA,alongB) {
  ## dimA[alongA] and dimB[alongB] must match.  All the submatrices corresponding to these
  ## dimensions are elementwise multiplied.  Important special cases:  with b a vector, diag(b) %*% A
  ## is the same as b * A , which in turn matches tensord(matrix(b, ncol=1), A, 1, 1) and A %*% diag(b) is
  ## is the same as t(b * t(A)), which in turn matches tensord(matrix(b,ncol=1), A, 1, 2).  Using the diag(b)
  ## expressions in these cases is very inefficient for large matrices.
  if (is.null(dim(A)) ) dim(A) <- length(A)
  if (is.null(dim(B))) dim(B) <- length(B)
  if (length(dim(B)) < length(dim(A)) {   #inefficient ordering.  Switch.
    scratch <- A
    A <- B
    B <- scratch
    scratch <- alongA
    alongA <- alongB
    alongB <- scratch
  }
  dimA <- dim(A)
  dimB <- dim(B)
  stopifnot(identical(dimA[alongA], dimB[alongB]))
  ndA <- length(dimA)
  ndB <- length(dimB)
  outA <- setdiff(1:ndA, alongA)
  outB <- setdiff(1:ndB, alongB)
  A <- aperm(A,c(alongA,outA))
  B <- aperm(B,c(alongB,outB))
  dim(A) <- c(prod(dimA[alongA]), prod(dimA[outA]))
  AB <- array(0, c(prod(dimA[alongA]), prod(dimA[outA]), prod(dimB[outB])))
  for (i in dim(AB)[2]) AB[ , i, ] <- A[ , i] * B
  dim(AB) <- c(dimA[alongA], dimA[outA], dimB[outB])
  return(AB)
}

  
