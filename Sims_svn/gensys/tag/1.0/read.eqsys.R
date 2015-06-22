read.eqsys <- function(file) {
  ix <- 0
  data <- readLines(con=file)
  ## remove comment lines
  data <- data[-grep("^[ \t]*##",data)]
  data <- data[-grep("^[ \t]*$", data)]
  ## so comment lines are any starting with 0 or more blanks and
  ## tabs, followed by ##.  Blank lines also ignored.
  neq <- length(data)/2 - 3
  name <- vector("character", neq)
  forward <- vector("logical", neq)
  eq <- vector("expression", neq)
  vlist <- vector("character")
  param <- vector("character")
  shock <- vector("character")
  data <- matrix(data, 2, neq + 3)
  for (iq in 1:neq) {
    cnct <- textConnection(data[1,iq])
    nf <- scan(cnct, what="char", quiet=TRUE)
    close(cnct)
    nc <- nchar(nf)
    forward[iq] <- identical(substr(nf,nc,nc), "*")
    if (forward[iq]) nf <- substr(nf, 1, nc - 1)
    ## asterisk character is a problem if we index by equation name
    name[iq] <- nf
    eq[iq] <- parse(text=data[2,iq])
  }
  cnct <- textConnection(data[2, neq+1])
  vlist <- scan(cnct, what="char", quiet=TRUE)
  close(cnct)
  cnct <- textConnection(data[2, neq+2])
  param <- scan(cnct, what="char", quiet=TRUE)
  close(cnct)
  cnct <- textConnection(data[2, neq+3])
  shock <- scan(cnct, what="char", quiet=TRUE)
  close(cnct)
  class(eq) <- c("eqsys", "expression")
  names(eq) <- name
  eq <- structure(eq, forward=forward, vlist=vlist, param=param, shock=shock)
  return(eq)
}
