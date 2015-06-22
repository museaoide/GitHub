\name{qzswitch}
\alias{qzswitch}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Interchanges adjacent eigenvalues in a qz decomposition }
\description{
The qz decomposition is in the form Q, A, B, Z, with A and B triangular
and Q and Z Hermitian (complex orthonormal).  qzswitch preserves this
structure while reordering an adjacent pair of diagonal rows and columns.
}
\usage{
qzswitch(i = 1, qzlist)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{i}{Location of first element to be interchanged }
  \item{qzlist}{List defining the transform as returned by \code{qz()}}
}
%\details{
%}
\value{
  ~Describe the value returned
  If it is a LIST, use
  \item{ }{Description of 'comp1'}
  \item{comp2 }{Description of 'comp2'}
  ...
}
\references{ ~put references to the literature/web site here ~ }
\author{ ~~who you are~~ }
\note{ ~~further notes~~ 

 ~Make other sections like Warning with \section{Warning }{....} ~
}
\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function(i=1,a=as.complex(1),b=as.complex(1),q=as.complex(1),z=as.complex(1))
  {
    ## Takes U.T. matrices a, b, orthonormal matrices q,z, interchanges 
    ## diagonal elements i and i+1 of both a and b, while maintaining 
    ## qaz' and qbz' unchanged.  If diagonal elements of a and b 
    ## are zero at matching positions, the returned a will have zeros at both 
    ## positions on the diagonal.  This is natural behavior if this routine is used 
    ## to drive all zeros on the diagonal of a to the lower right, but in this case 
    ## the qz transformation is not unique and it is not possible simply to switch 
    ## the positions of the diagonal elements of both a and b. 
    realsmall <- 1e-7; 
    ##realsmall<-1e-3;
    q <- t(Conj(q))                     #This is needed because the code was originally for matlab, where it is q'az' that
                                        # is preserved.
    A <- a[i,i]; d <- b[i,i]; B <- a[i,(i+1)]; e <- b[i,(i+1)]
    c <- a[i+1,i+1]; f <- b[i+1,i+1]
    ## a[i:(i+1),i:(i+1)]<-[A B; 0 c]
    ## b[i:(i+1),i:(i+1)]<-[d e; 0 f]
    if (abs(c)<realsmall & abs(f)<realsmall)
      {
        if (abs(A)<realsmall)
          {
            ## l.r. coincident 0's with u.l. of a<-0; do nothing
            return(list(a=a,b=b,q=q,z=z))
          } else
        {
          ## l.r. coincident zeros; put 0 in u.l. of a
          wz <- c(B, -A)
          wz <- wz/sqrt(sum(conj(wz)*wz))
          wz <- array(c(wz,Conj(wz[2]),-Conj(wz[1])),dim=c(2,2))
          xy <- diag(2)
        }
      } else
    {
      if (abs(A)<realsmall && abs(d)<realsmall)
        {
          if (abs(c)<realsmall)
            {
              ## u.l. coincident zeros with l.r. of a<-0; do nothing
              return(list(a=a,b=b,q=q,z=z))
            } else
          {
            ## u.l. coincident zeros; put 0 in l.r. of a
            wz <- diag(2)
            xy <- c(c,-B)
            xy <- xy/sqrt(sum(xy*Conj(xy)))
            xy <- t(matrix(c(Conj(xy[2]), -Conj(xy[1]),xy),nrow=2,ncol=2))
          }
        } else
      {
        ## usual case
        wz <- c(c*e-f*B, Conj(c*d-f*A))
        xy <- c(Conj(B*d-e*A), Conj(c*d-f*A))
        n <- sqrt(wz \%*\% Conj(wz))
        m <- sqrt(xy \%*\% Conj(xy));
        if (Re(m)<1e-12*100)
          {
            ## all elements of a and b proportional
            return(list(a=a,b=b,q=q,z=z))
          }
        wz <- wz/n
        xy <- xy/m
        wz <- matrix(c(wz, -Conj(wz[2]),Conj(wz[1])),byrow=TRUE,ncol=2,nrow=2)
        xy <- matrix(c(xy,-Conj(xy[2]), Conj(xy[1])),byrow=TRUE,ncol=2,nrow=2)
      }
    }
    a[i:(i+1),] <- xy \%*\% a[i:(i+1),]
    b[i:(i+1),] <- xy \%*\% b[i:(i+1),]
    a[,i:(i+1)] <- a[,i:(i+1)] \%*\% wz
    b[,i:(i+1)] <- b[,i:(i+1)] \%*\% wz
    z[,i:(i+1)] <- z[,i:(i+1)] \%*\% wz
    q[i:(i+1),] <- xy \%*\% q[i:(i+1),]
    q <- t(Conj(q))
    return(list(a=a,b=b,q=q,z=z))
  }
}
\keyword{ ~kwd1 }% at least one, from doc/KEYWORDS
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
