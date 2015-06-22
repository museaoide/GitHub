discreteCAPM <- function(x,lambda){
  ### x is a matrix that determines the joint pdf of R and s, The objective is to max
  ### R*s-.5*(R*s)^2 - lambda*(mutual info between R and s)
