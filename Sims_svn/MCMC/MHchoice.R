MHchoice <- function(lhprop, lhtarget, theta0, theta1) {
  rho <- lhtarget(theta1) * lhprop(theta0)/(lhtarget(theta0) * lhprop(theta1))
  if( runuf(1) > rho ) return (theta1) else return(theta0)
}
