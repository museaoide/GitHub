invTime <- function(tsdate, tsobj) {
  ## for tsdate entered as c(1983,3) or as 1983.5, returns the index corresponding to that
  ## date in tsobj.  tsdate can be an n x 2 array (in which case it contains c(1983,3) type
  ## dates, or it can be a vector.  A vector with two elements can create ambiguity, but the
  ## it is usually resolvable. To be safe, use 1 x 2 matrices
  ## Note that this is checking equality of Floating point numbers, so it relies on real number
  ## arguments being computed exactly as does time(), e.g. 1983+2/12, not 1972.167.
  ##--------------------------
  start <- tsp(tsobj)[1]
  end <- tsp(tsobj)[2]
  freq <- tsp(tsobj)[3]
  if (is.null(dim(tsdate)) && length(tsdate)==2) { # ambiguity
    d2 <- tsdate[1] + (tsdate[2] - 1) / freq
    d2OK <- d2 >= start   && d2 <= end 
    d1OK <- all((tsdate >= start ) & (tsdate <= end) )
    if (d2OK && d1OK) {
      warning("ambiguous date")
      return(match(d2, time(tsobj)))
    } else {
      ifelse( d2OK, match(d2, time(tsobj)), match(tsdate, time(tsobj)))
    }
  } else {
    if (is.null(dim(tsdate))) {
      return(match(tsdate, time(tsobj)))
    } else {
      tsd <- tsdate[ ,1] + (tsdate[ ,2] -1)/freq
      return(match(tsd, time(tsobj)))
    }
  }
}
