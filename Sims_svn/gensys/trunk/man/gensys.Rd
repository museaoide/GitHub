\name{gensys}
\alias{gensys}
\title{ Solve linear rational expectations models }
\description{
Solve a first order linear RE model using the QZ (i.e. gneralized Schur)
decomposition.  Also test for existence and uniqueness of solutions.
The model is in the form \deqn{\Gamma_0 y_t = \Gamma_1 y_{t-1} + c_0 + \Psi z_t + \Pi \eta_t}
with \eqn{E_t[\eta_{t+1}]=0}.
}
\usage{
gensys(g0, g1, c0, psi, pi, div = -1)
}
\arguments{
  \item{g0}{ contemporaneous coefficient matrix \eqn{\Gamma_0}}
  \item{g1}{ lagged coefficient matrix \eqn{\Gamma_1}}
  \item{c0}{ constant vector}
  \item{psi}{ coefficient \eqn{Psi} on exogenous \eqn{z_t} }
  \item{pi}{ coefficient \eqn{Pi} on expectational errors \eqn{\eta_t} }
  \item{div}{ roots exceeding this in absolute value are not allowed to appear in the solution }
}
%\details{
%  ~~ If necessary, more details than the description above ~~
%}
\value{
	Solution is in the form 
	\deqn{y_t = \theta_c + \theta_1 y_{t-1} + \theta_0 z_t +\sum_{s=1}^\infty \theta_y \theta_f^s \theta_z z_{t+s}}
  \item{G1}{\eqn{\theta_1}}
  \item{C}{\eqn{\theta_c}}
  \item{impact}{\eqn{\theta_0}}
  \item{fmat}{\eqn{\theta_f}}
  \item{fwt}{\eqn{\theta_z}}
  \item{ywt}{\eqn{\theta_y}}
  \item{gev}{generalized eigenvalues}
  \item{eu}{\code{eu[1]} = solution exists; \code{eu[2]} = solution is unique;
  For the \code{eu} values, 0 means \code{FALSE}, 1 means \code{TRUE}.  (0,0) is possible.  It means
  there is no stable solution, and even if a stable solution were imposed by forcing constraints
  on the behavior of \eqn{z}, there would be indeterminacy.  (-2,-2) is possible.  It means there
  are coincident zeros in \code{gev}, which implies that the system is not of full rank --- some equations
  are linear combinations of others.}
  \item{loose}{When \code{eu[2]==0}, the columns of \code{loose} span the space of indeterminacies.  Arbitrary
  serially uncorrelated shocks to the equation system in these directions cannot be ruled out.
  }
}
\author{ Chris Sims <sims@princeton.edu> }
%\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\section{Example model}{
The example code solves this model:\cr
\tabular{ll}{
Fisher equation: \tab \eqn{i_t =  E_t p_{t+1} - p_t}\cr
policy: \tab    \eqn{i_t = \rho i_{t-1} + \theta p_{t-1} + \epsilon_t}
}
}
\examples{
theta <- .1
rho <- .9
g0 <- matrix(c(0,-1, 1,0), 2, 2, byrow=TRUE)
g1 <- matrix(c(-1,-1, rho, theta), 2, 2, byrow=TRUE)
Pi <- c(1,0)
Psi <- c(0,1)
c0 <- c(0,0)
gout <- gensys(g0, g1, c0, Psi, Pi)
gout
rm(gout)
}
\keyword{math} 
