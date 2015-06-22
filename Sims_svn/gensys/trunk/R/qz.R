qz <- function(a=array(1,dim=1,1),b=array(1,dim=1,1)) {
  ## R does not load all of lapack in its own lapack.so library.  The routine
  ## needed here, zgges, therefore has to be loaded directly.
  if(!is.loaded("zgges")){
    filename <- "/usr/lib/liblapack.so" #not working with atlas-compiled R.
    if (file.access(filename)) 
      filename <- "/usr/lib/liblapackgf-3.so"
    dyn.load(filename, now=FALSE)
  }
  N<-dim(a)[1]
  SDIM<-as.integer(1);
  ALPHA<-vector("complex",N)
  BETA <-vector("complex",N)
  q <- array(data=0+0i,dim=c(N,N))
  z <- array(data=0+0i,dim=c(N,N))
  LWORK <- as.integer(33*N)
  WORK <- vector("complex",length=LWORK)
  RWORK <- vector("numeric",length=8*N)
  INFO <- 0
## -------------------------------
## Notes:  The reordering could be done inside zgges.f.  This would require writing a fortran module that contained both the
## routine to check the size or sign of the real part of the root and a separate program, that sets div in the module.  (All this to
## get arund the fact that zgges wants the comparison function to have just two arguments.)  Also, to implement this in f77, and thereby
## guarantee wider portability, one would have to use a common block and ENTRY, rather than a module and double,save.  
  ##-------------------------------
  ## Returns a,b,q,z s.t. q %*% t(Conj(q)) = z %*% t(Conj(z)) = I, a and b upper triangular,
  ## q %*% a %*% t(Conj(z)) = A, and q %*% b %*% t(Conj(z))=B.  The diagonal of b should be
  ## non-negative and real (a normalization).  gev is a two-column array containing the generalized
  ## eigenvalues of the pair A,B.
  ##
  ## rc=0 indicates return without problems.
  ## rc=-18 indicates not enough work space
  ## was allocated.  The current code attempts to be overgenerous in allocating workspace
  ## based on empirical testing of the lapack code with random matrices.  If rc=-18 ever now
  ## occurs in practice, the code could be modified to first test for required space.
  ## Also, if the current version's space requirements are a problem, the code can surely
  ## be modified to use less workspace, with some possible loss in efficiency.
  ## rc=-j indicates an illegal argument type for the j'th argument to the lapack routine.  This should
  ## not occur.
  ## -------------------------------------
  ## The call to dyn.load works on  SuSE 9 and 10 systems with the distribution's installation
  ## of lapack.  If your lapack is at a different location, modify accordingly.
  ## -------------------------------------
  ## The first two V's indicate that we want q and z computed.  The quoted
  ## "N" indicates we want no sorting.  (This would require producing a
  ## different fortran routine for each possible choice of the "stake" around
  ## which we are sorting.) The first dum is where the sort criterion routine's
  ## name would go.  N is the 
  ## dimension of all the matrices.  SDIM would be the number of eigenvalues
  ## for which the sort criterion is true, but since we're not sorting it
  ## will come back 0, always.  ALPHA and BETA are the generalized eigenvalue
  ## vectors (the diagonals of the returned a and b).  WORK and RWORK
  ## are workspace, the former complex, the latter real.  The last dum is
  ## a placeholder for logical workspace needed only if a sort is done.
  ## INFO is 0 if all goes well, -j if the j'th argument has an incorrect
  ## form, j>0 if the algorithm failed, but ALPHA(k) and BETA(k) are accurate
  ## for k>j, N+1 if there was another sort of failure.
  ## ------------------------------
  out<-.Fortran("zgges","V","V","N","dum",N,as.complex(a),N,as.complex(b),N,
                SDIM,ALPHA,BETA,q,N,z,N,WORK,LWORK,RWORK,"dum",INFO)
  gev<-matrix(c(out[[11]],out[[12]]),nrow=N,ncol=2)
  ## if you run into problems with workspace, uncomment the line below and use the
  ## message it produces to modify the setting of LWORK above.
  ## cat("workspace needed:",out[[17]][1],"\n")
  return(list(a=matrix(out[[6]],nrow=N,ncol=N),
              b=matrix(out[[8]],nrow=N,ncol=N),
              q=matrix(out[[13]],nrow=N,ncol=N),
              z=matrix(out[[15]],nrow=N,ncol=N),
              gev=gev,
              rc=out[[21]]))
}
     
     
     
     
  
