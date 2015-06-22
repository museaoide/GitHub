\name{qzdiv}
\alias{qzdiv}
\title{Sort roots of QZ transform into large and small  }
\description{
  Rearranges rows and columns of the output of \code{qz} so that large
  roots are at the bottom.
}
\usage{
qzdiv(stake, qzlist)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{stake}{roots larger than this are movd to the bottom }
  \item{qzlist}{output list from \code{qz}; the \code{rc} entry in
    the list is not used or modified }
}

\value{
The function returns a list whose elements \code{a,b,q,z,gev} are the
re-ordered input matrices and generalized eigenvalues.  The returned
list is of the same form as the input list.
}
\author{Christopher Sims, \email{sims@princeton.edu}, \url{http://www.princeton.edu/~sims} }
\examples{
A <- matrix(c(1,1,1, 1,2,1, 0,1,3),3,3)
B <- matrix(c(3,0,0, 1,2,1, 1,1,1),3,3)
A
B
qzlist <- qz(A,B)
qzlist
qzlist$q %*% qzlist$a %*% Conj(t(qzlist$z))
qadiv(1,qzlist)
qzlist
qzlist$q %*% qzlist$a %*% Conj(t(qzlist$z))
rm(A,B,qzlist)
}
\keyword{math }% at least one, from doc/KEYWORDS
