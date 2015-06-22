\name{gensys-package}
\alias{gensys-package}
\alias{gensys}
\docType{package}
\title{
Solve Linear Rational Expectations Models
}
\description{Solves linear rational expectations models and provides lower level
tools that can be useful in analyzing the structure of a difficult model.
}
\details{
\tabular{ll}{
Package: \tab gensys\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2006-05-06\cr
License: \tab LGPL\cr
}
\code{gensys} is the model solver.  \code{qz} is an \R wrapper for a \code{lapack}
QZ decomposition routine.  The raw QZ decomposition can be useful for diagnosing
problems with a model.  \code{qzdiv} sorts the generalized eignvalues into larger and
smaller sets.  \code{qzswitch} is used by \code{qzdiv} and is probably seldom useful by
itself.  \code{gstate} is an aid to interpreting the output of \code{gensys}, attempting to
make dimension-reducing simplifications by identifying candidate state
vectors.  \code{gensysct,qzdivct} are for continuous time models.
}
\author{
Chris Sims: \email{sims@princeton.edu}, \url{http://www.princeton.edu/~sims}

}
\references{
Sims, Christopher A., \dQuote{Solving linear rational expectations models}.
\emph{Computational Economics}, 20(1-2):1-20, 2001.
}
\keyword{ package }
