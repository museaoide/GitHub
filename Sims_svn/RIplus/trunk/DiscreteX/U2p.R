U2p <- function(x,y) {
    a <- .5                             #a is cross-effect in demand function
    prf <- (1 - x[1] + a * x[2]) * (x[1] - y[1]) + (1 - x[2] + a * x[1]) * (x[2] - y[2])
    return(prf)
}
