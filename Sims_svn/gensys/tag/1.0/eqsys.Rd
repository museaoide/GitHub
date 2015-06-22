\name{eqsys}
\alias{eqsys}
\alias{print.eqsys}
\alias{read.eqsys}
\alias{c.eqsys}
\title{eqsys S3 class}
\description{A system of equations object}
\usage{print.eqsys(eq)
\method{read}{eqsys}(con=stdin())
\method{c}{eqsys}(\dots, recursive=FALSE)
}
\arguments{
  \item{eq}{an \code{eqsys} object.}
  \item{con}{a connection, from which the system is read.}
  \item{\dots}{a comma-separated set of \code{eqsys} objects to be concatenated.}
}
\value{
  \item{print}{a text representation of the system suitable for
    interactive printout.  Forward-looking equations marked with an asterisk.}
  \item{read}{an \code{eqsys} object}
  \item{c}{an eqsys object that stacks up the arguments}
}
\details{An \code{eqsys} object has equations stored as a vector of
  expressions that should evaluate to zero (no equal signs).  It also
  has attributes \code{names}, a character vector that names the
  equations, \code{forward}, a logical vector that is \code{TRUE} in
  the position corresponding to each forward-looking equation, \code{vlist}, a character vector of
 variable names, \code{param}, a character vector of parameter names, and
 \code{shock}, a character vector of shock names.  A dynamic system will
 generally also have lagged variables or derivatives.  These should be
 entered with \code{l} or \code{dot} suffixes (e.g. \code{gdpl} or
 \code{gdpdot}), and are not recorded separately in the \code{eqsys}
 object.
 
 The input for \code{read} should start, if there are n equations, with
   2n line pairs.   The first line in each pair contains the equation name, followed
   by an asterisk if th equation is forward-looking.  The second line contains the expression defining the
   equation.  After these 2n lines should come a single line containing
   the string 'vlist', followed by a single line with the variable
   names, separated by blanks.  Then there should be a single line
   containing 'param', followed by a single line parameter names,
   separated by blanks, followed by a single line containing 'shock',
   followed by a single line containing the shock names.

   Bugs:  Should allow multiline equation expressions, vlists, etc., with blank line separators for equaitons.  Should (?) specify
   which variables are pegged and equations dropped in steady state calculation (because of unit roots).  
}
\references{
}
\seealso{
}
\examples{
}
\author{}
\keyword{}
