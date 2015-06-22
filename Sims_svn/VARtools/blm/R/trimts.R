#' Trim a multiple time series object to eliminate NA's
#' 
#' Using \code{cbind} to combine lagged values of time series
#' creates NA's at the start and end.  This can be avoided by
#' using \code{ts.intersect} instead of cbind.  But if one starts
#' with an mts object with series of different time spans, the
#' NA's will be there, and ts.intersect ignores them.  
#' 
#' @param y The mts object to be trimmed.
trimts <- function(y) {
  if(!is.ts(y)) error("arg to trimts() not a time series")
  if(is.null(dim(y))) dim(y) <- c(length(y), 1)
  firstGood <- match(TRUE, apply(y, MARGIN=1, function(x) !any(is.na(x))))
  firstGood <- time(y)[firstGood]
  y <- window(y, start=firstGood)
  lastGood <- match(TRUE, apply(y, MARGIN=1, function(x) any(is.na(x))))
  if(!is.na(lastGood)) {
    lastGood <- time(y) [lastGood-1]
    y <- window(y, end=lastGood)
  }
  return(y)
}