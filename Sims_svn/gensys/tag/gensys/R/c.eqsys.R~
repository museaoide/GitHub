c.eqsys <- function(...,recursive=FALSE) {
  ## still need to remove duplicates in vlist, param, shock
  args <- list(...)
  forward <- vector("logical")
  vlist <- vector("character")
  param <- vector("character")
  shock <- vector("character")
  for (eq in args) {
    forward <- c(forward, attr(eq, "forward"))
    vlist <- union(vlist, attr(eq, "vlist"))
    param <- union(param, attr(eq, "param"))
    shock <- union(shock, attr(eq, "shock"))
  }
  structure(NextMethod(...), class=c("eqsys","expression"), forward=forward, vlist=vlist,
            param=param, shock=shock)
}
