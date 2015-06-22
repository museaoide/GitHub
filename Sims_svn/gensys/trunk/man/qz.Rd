\name{qz}
\alias{qz}
\alias{qzswitch}
\alias{qzdiv}
\alias{qzdivct}
\title{ QZ (generalized Schur) decomposition }
\description{
 \code{qz} Calculates the QZ decomposition of a pair of square
 matrices.  \code{qzdiv} sorts the eigenvalues of the decomposition,
 based on the absolute value of the eigenvalues.  \code{qzdivct} sorts
 based on the real part of the eigenvalues.  \code{qzswitch} is used by
 the two previous programs; it interchanges two adjacent eigenvalues.
}
\usage{
qz(a = array(1, dim = 1, 1), b = array(1, dim = 1, 1), div = 1 + 1e-14)
qzdiv(stake, qzlist)
qzdivct(stake, qzlist)
qzswitch(i, qzlist)
}
\arguments{
  \item{a}{ left hand side input matrix }
  \item{b}{ right hand side input matrix }
  \item{stake}{Sort divider.  \code{qzdiv} puts roots larger than this
    in absolute value in upper left.  \code{qzdivct} puts roots with
    real part larger than this in upper left }
  \item{qzlist}{a list of the components of the qz decomposition}
  Components of \code{qzlist}:
  \item{a}{left triangular matrix}
  \item{b}{right triangular matrix}
  \item{q}{left unitary matrix}
  \item{z}{right unitary matrix}
}
\details{  If \eqn{A} and \eqn{B} are two \eqn{N x N}  matrices, we can find
  \eqn{a,b,q,z} with \eqn{A = q a z'}, \eqn{B = q b z'}, \eqn{a,b} both
  upper triangular, and \eqn{q'q = z'z = I}, where in these formulas the
  \sQuote{'} indicates both transposition and complex conjugation.  The
  computation is done by a call to the lapack routine \code{zgges}.  [The
  \code{qzdiv} sorting could be done more efficiently within
  \code{zgges}, but this would require writing fortran code to implement
  the sort comparisons.]  
}

\value{
  \item{a,b,q,z}{the components of the qz decomposition as described above.}
\item{gev}{ \code{N x 2} matrix of generalized eigenvalue pairs.}
  \item{rc}{return code}
   The return code values are
  \item{0:}{all went well}
  \item{\code{j > 0}:}{the algorithm failed, but
    \code{gev[k,]} is accurate for \code{k > j}} 
  \item{\eqn{N+1}:}{some other sort of failure (\eqn{N} is the
    dimension of the matrices)}
  \item{-18:}{workspace overflow; see comments in the source if this
    occurs.}
}

\references{Golub, Gene H. and Charles F. Van Loan, \emph{Matrix
    Computations}, 3d edition, Johns Hopkins University Press 1996,
  section 7.7}
\author{Christopher Sims, \email{sims@princeton.edu}, \url{http://www.princeton.edu/~sims }}
\examples{
qzout <- qz(matrix(c(1,0,0, 1,0,0, 1,1,1),3,3),matrix(c(0,1,0, 1,1,0,
    1,1,0),3,3))
qzout
qzdiv(stake=1.1, qzout)
rm(qzout)
}
\keyword{math }% at least one, from doc/KEYWORDS
