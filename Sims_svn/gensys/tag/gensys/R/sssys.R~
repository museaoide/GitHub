#' sssys
#'
#' Transform an \code{eqsys} dyanmic system to a static system for use
#' in finding the steady state, usually with \code{ssSolve()}.
#' 
#' @param eq an eqsys object with a declared vlist attribute and with
#'   lags denoted by l's attached to the end of the variable names.
#'   (Notice that this means that if a variable name ends in l, it
#'   should not be another variable name when the terminal l is removed.)
#' @return A new system with all lagged variable names
#'   replaced by current variable names.
#' 
sssys <- function(eq) {
  nq <- length(eq)
  nv <- length(attr(eq, "vlist"))
  eq2 <- eq
  for ( iv in 1:nv) {
    eq2 <- parse(text=gsub(paste("(\\b", attr(eq,"vlist")[iv], ")l\\b", sep=""),  "\\1", eq2, perl=TRUE))
  }
  attributes(eq2) <- attributes(eq)
  return(eq2)
}
