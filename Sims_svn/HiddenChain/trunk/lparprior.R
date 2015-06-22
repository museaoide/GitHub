lparprior <- function(parvec, y0, parPrior) {
  ## parvec has components corresponding to parameters that don't change with state, and also
  ## components that do change.  In this version we consider a reduced form VAR with time-varying
  ## variance matrix, modeled as Sigma = sum_j sij^2 * vj %*% t(vj), where the variance weights on
  ## the components vary with the state, while the component vectors vj do not.  The vj's have
  ## to be normalized to have length one.  
  ##--------------------- make prior proper------------------------
  ny <- dim(y0)[2]
  nlag <- dim(y0)[1]
  with(parvec,
       with(parPrior,{
         
}

