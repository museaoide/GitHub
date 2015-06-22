#' eqsys
#'
#' A system of equations object
#'
#' @param eq an \code{eqsys} object
#' @param \item{file}{a connection (e.g. a file name), from which the system 
#'  is read.}
#'  @param \dots a comma-separated set of \code{eqsys} objects to be
#'     concatenated.
#'
#'@return
#' \list{
#'  \item{print}{a text representation of the system suitable for
#'    interactive printout, but also usable as input to \code{read.eqsys()}. }
#'  \item{read}{an \code{eqsys} object}
#'  \item{c}{an eqsys object that stacks up the arguments}
}
#' @details
#' An \code{eqsys} object has equations stored as a vector of  expressions
#' that should evaluate to zero (no equal signs).  It also
#' has attributes \code{names}, a character vector that names the
#' equations, \code{forward}, a logical vector that is \code{TRUE} in
#' the position corresponding to each forward-looking equation, \code{vlist},
#' a character vector of variable names, \code{param}, a character vector
#' of parameter names, and  \code{shock}, a character vector of shock names.
#' A dynamic system will  generally also have lagged variables or derivatives.
#' These should be  entered with \code{l} or \code{dot} suffixes
#' (e.g. \code{gdpl} or  \code{gdpdot}), and are not recorded separately in
#' the \code{eqsys}  object.
#' 
#' The input for \code{read} should start, if there are n equations, with 2n
#' line pairs.   The first line in each pair contains the equation name,
#' followed by an asterisk if the equation is forward-looking.  The second
#' line contains the expression defining the equation.  After these 2n lines
#' should come a single line containing the string 'vlist', followed by a
#' single line with the variable names, separated by blanks.  Then there
#' should be a single line containing 'param', followed by a single line of
#' parameter names, separated by blanks, followed by a single line containing
#' 'shock', followed by a single line containing the shock names.  For a model
#' with no shocks, the shock name line should contain "NONE". Blank lines
#' are ignored, and any line beginning with ## or blanks and tabs followed
#' by ## is ignored (i.e. is a comment line).
#'
#'   Bugs:  Should allow multiline equation expressions, vlists, etc., with
#' blank line separators for equaitons.  Should (?) specify which variables
#' are pegged and equations dropped in steady state calculation (because of
#' unit roots).  
#'
#' @references Sims, C. A., "Solving Linear Rational Expectations Models"
#' \itshape{Computational Economics}, 2002, 20, p, 1-20
#'
#' @seealso{\code{gensys}
#' @examples
#' ## Simple two-equation asset pricing
#' ## equation system specification file to be read by read.eqsys:  (without
#' initial ##'s)
#' ## -------------------------
#' ## divDynamics
#' ## div - (.1 + .9 * divl + eps)
#' ##
#' ## pricing*
#' ## pl - beta * (p + div)
#' ##------------------------
#' ## vlist
#' ## div p
#' ##
#' ## param
#' ## beta
#' ##
#' ## shock
#' ## eps
#' 
#' @describeIn eqsys Print an \code{eqsys} object
print.eqsys <- function(eq)
#' @describeIn eqsys read a text file containing an \code{eqysys} object
read.eqsys <- function(file)
#' @describeIn eqsys concatenate \code{eqsys} objects
c.eqsys <- function(\dots, recursive=FALSE)
