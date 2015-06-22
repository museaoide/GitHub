U1.1 <- function (a, b) 
{
    .expr1 <- a - b
    .expr2 <- .expr1^2
    .value <- -.expr2^0.55
    .grad <- array(0, c(length(.value), 1L), list(NULL, c("a")))
    .grad[, "a"] <- -(0.55 * (2 * .expr1 * .expr2^-0.45))
    .grad[a==b, "a"] <- 0
    .value <- ifelse(a==b, 0, .value)
    attr(.value, "gradient") <- .grad
    .value
}
