g0g1evalc <- function(dexpr, xd=vector("numeric",length(x)), x, shock, experr, param) {
  ##  dexpr:     derivative expression vector, as produced by g0g1dc
  ##      x:     nx by 1 matrix of values of x, names in dimnames(x)[[1]]
  ##     xd:     nx by 1 matrix of values of xdot, names in dimnames(x)[[1]]
  ##  shock:     exogenous shocks (all assumed centered at 0)
  ## experr:     vector of numbers or names of equations with expectational errors
  ##  param:     vector of values of parameters (all other symbols) in the expressions
  ## NB x, xd, and param must all be vectors with named components, or lists.
  if (is.null(names(x))) { names(x) <- dimnames(x)[[1]]}
  if (is.null(names(xd))) names(xd) <- paste(names(x),"dot",sep="")
  if (is.null(names(shock))) names(shock) <- dimnames(shock)[[1]]
  if (is.null(names(param))) names(param) <- dimnames(param)[[1]]
  nex <- length(dexpr)
  nx <- length(x)
  nxd <- length(xd)
  nshock <- length(shock)
  g0 <- matrix(0, nex, nxd )
  g1 <- matrix(0, nex, nx )
  Psi <- matrix(0, nex, nshock)
  Pi <- matrix(0, nex, length(experr))
  for ( i in 1:nex ) {
    ## browser()
    dvec <- attr(eval(eval(dexpr[i]), as.list(c(xd,x,shock,param))),"gradient")
    loc <- 0
    g0[i,] <- dvec[1:nx]
    loc <- nx
    g1[i,] <- dvec[(loc+1):(loc+nxd)]
    loc <- loc + nxd
    Psi[i,] <- dvec[(loc+1):(loc+nshock)]
  }
  Pi[experr,] <- diag(length(experr))
  dimnames(g0) <- list(names(dexpr),names(x))
  dimnames(g1) <- dimnames(g0)
  if(is.numeric(experr)) experr <- names(dexpr)[experr]
  dimnames(Psi) <- list(names(dexpr),names(shock))
  dimnames(Pi) <- list(names(dexpr), experr)
  return(list(g0=g0, g1=-g1, Psi=-Psi, Pi=-Pi))
}
