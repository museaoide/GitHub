#' Step responses from ts regression coefficients
#'
#' Responses to persistent unit perturbations in right-hand-side variables
#' 
#' Applies \code{filter} twice, using first coefficients as lagged dependent variable
#' coefficients in a recursive filter and the coefficients on variable nv in a 
#' one-side convolution filter.
#'  
#' @param bdraw matrix with each row a draw of the coefficient vector
#' @param lags a vector, with each element the number of lags of the corresponding
#'        variable that appears in in the coefficient vector
#' @param to rhs variable number for which we want the response
#' @param horiz number of periods in the calculated responses
#' @return matrix with each column a draw of the response
tsregImpDraw <- function(bdraw, lags, to, horiz) {
  dvlag <- lags[1]
  tolag <- lags[to]
  ndraw <- dim(bdraw)[1]
  cfldv <- bdraw[ , 1:dvlag]
  cfto <- bdraw[ , sum(lags[1:to]) + 1:tolag]
  resp <- apply(X=cfto,  MARGIN=1, FUN=filter, x=c(rep(0,tolag), rep(1,horiz)), method="conv", sides=1)
  for (idraw in 1:ndraw) {
    resp[-(1:tolag) , idraw] <- filter(resp[-(1:tolag), idraw], filter=cfldv[idraw, ], method="recurs")
  }
  return(resp[-(1:tolag), ])
}