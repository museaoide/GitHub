#' Bayesian linear regression
#' 
#' Regression with a proper conjugate prior specified with dummy observations.
#' 
#' With dummy observations and degrees of freedom and scale for
#' the marginal distribution of disturbance variance, constructs
#' estimates and provides the marginal data density, allowing
#' model comparison via Bayes factors.
#' @param x The observed right-hand side variable data
#' @param y The observed left-hand side variable data
#' @param xp The dummy observation matrix for the right-hand side
#' @param yp The dummy observations for the left-hand side
#' @param dfp Degrees of freedom for the inverse-gamma prior on residual 
#' variance
#' @param scalepv Scale for prior on variance
#' @return List with elements 
#' \itemize{
#'   \item lsout:  The returned output from \code{\link{lm.fit}}..
#'   \item postdf: Degrees of freedom in the inverse-gamma posterior for 
#'         residual variance.
#'   \item postssr: Sum of squared residuals, over real and dummy observations,
#'         at posterior mode/mean. 
#'   \item lmdd: log of marginal data density.
#'   \item vcv: Scale matrix for posterior \eqn{t} distribution on parameters.
#'              Variance matrix if multiplied by \code{postdf/(postdf-2)}
#'   \item sdcoeff: Vector of standard deviations of coefficients
#'         \code{sqrt(diag(vcv))}
#' }
#' @details Constant term must be introduced explicitly as a column of ones in 
#'    \code{x}.  With \code{dfp=1}, the inverse of residual variance has mean
#'    \code{1/scalevp}.  \code{scalevp} must be
#'    increased in proportion to \code{dfp} to keep the center of the distribution
#'    for residual variance stable as it becomes more precise. 
#' @seealso \code{\link{tsregPrior}}, which prepares a symmetric prior, biased toward
#'    persistence, for models with lags.  \code{\link{lagts}} for preparing 
#'    \code{mts} objects that can play the role of \code{x} in this program.
blm <- function(x, y, xp, yp, dfp=1, scalepv) {
  if (is.null(dim(y)))
  {
    y <- matrix(y, ncol=1)
  }
  if (is.null(dim(yp))) {
    yp <- matrix(yp, ncol=1)
  }
  stopifnot(dim(x)[1] == dim(y)[1], dim(xp)[1] == dim(yp)[1], dim(x)[2] == dim(xp)[2])
  X <- rbind(x, xp)
  Y <- rbind(y, yp)
  N <- dim(X)[1]
  k <- dim(X)[2]
  lsout <- lm.fit(X, Y)
  postdf <- (N-k)/2 + dfp
  postssr <- sum(lsout$residuals^2)
  lmdd <- lgamma(postdf) - lgamma(dfp) + dfp * log(scalepv) - postdf * log(.5 * postssr + scalepv)
                                - (postdf - dfp) * log(2 * pi) + .5 * determinant(crossprod(xp))$modulus
                                - .5 * determinant(crossprod(X))$modulus
  vcv <- solve(crossprod(qr.R(lsout$qr))) * postssr/(2 * postdf)
  sdcoeff <- sqrt(diag(vcv))
  return(list(lsout=lsout, postdf=postdf, postssr=postssr, lmdd=lmdd, vcv=vcv, sdcoeff=sdcoeff))
}
