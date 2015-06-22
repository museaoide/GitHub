HCsmooth <- function(psGYfilt=matrix(0,5,2), pyGs=function(y=0, S=c(1, 2), ...){}, transMat=diag(2),ydata,xdata,...) {
  ## psGYfilt:    filtered state probabilities, a nT x nState matrix, p(S_t | Y_t)
  ## pyGs:        function delivering pdf of y(t) given past y and current S_t
  ## transMat:    Markov transition matrix.  Columns add to one.
  ## ydata:       nT x ny data matrix
  ## xdata:       nT x nx exogenous data matrix (including column of 1's, if there are constant terms).
  ## lags:        number of y lags
  ##
  ## Note that The filtered state probabilities include an initial probability vector for it=lags.  This, too, gets updated here.
  psGYsmooth <- psGYfilt
  nT <- dim(psGYfilt)[1]
  nState <- dim(psGYfilt)[2]
  for (it in seq(nT-1,lags,-1)) {
    psGYsmooth[it, ] <- psGYfilt[it, ] * (( (psGYsmooth[it+1, ] / psGYfilt[it+1, ]) *
              exp(pyGs( ydata[it+1], 1:nState, yl=ydata[it:(it-lags+1)], x=xdata[it+1, ], ... )) )
                                          %*% transMat)
  }
  psGYsmooth <- psGYsmooth / as.vector(psGYsmooth %*% c(1,1,1))
  return(psGYsmooth)
}
