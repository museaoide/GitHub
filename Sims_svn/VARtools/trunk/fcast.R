fcast <- function(y0, By, Bx, xdata=NULL, const=TRUE, horiz, shocks=NULL) {
### By: equations x variables x lags
### Bx: equations x nx
### xdata: (lags+horiz) x nx
### y0: lags x nvar.  Must be a matrix, even if univariate.
### produces horiz-step-ahead forecast.
### Value is:
### yhat: (horiz+lags) x nvar
  if (is.null(dim(y0)))
    lags <- length(y0)
  else
    lags <- dim(y0)[1]
  if( is.null(shocks)) shocks <- matrix(0, horiz, dim(y0)[2])
  stopifnot(all.equal(dim(shocks), c(horiz,dim(y0)[2])))
  stopifnot( lags == dim(By)[3] )
  stopifnot( is.null(xdata) || (is.null(dim(xdata)) && length(xdata) == horiz+lags) || horiz+lags == dim(xdata)[1] )
  if (const) {
    if (is.null(xdata))
      xdata <- matrix(1,horiz+lags,1)
    else
      xdata <- cbind(xdata, matrix(1, horiz+lags, 1))
  } else {                              #no constant, or it's explicitly in xdata
    if (!is.null(xdata)) {
      if (is.null(dim(xdata))) xdata <- matrix(xdata,ncol=1)
    } else {                            #no constant, no xdata.  0 x 1  as placeholder
      xdata <- matrix(0, horiz+lags, 1)
      Bx <- matrix(1, 1, 1)
    }
  }
  if (is.null(dim(y0))) dim(y0) <- c(length(y0),1)
  nvar <- dim(y0)[2]
  nx <- dim(Bx)[2]
  yhat <- matrix(0,horiz+lags,nvar)
  yhat[1:lags,] <- y0
  Bmat <- aperm(By,c(3,2,1))            #lags by vbls by eqns
  Bmat <- Bmat[seq(lags,1,by=-1),,]     #reverse time index
  dim(Bmat) <- c(lags*nvar,nvar)
  for (it in 1:horiz){
    ydata <- yhat[it:(it+lags-1),,drop=FALSE]
    yhat[lags+it,] <- apply(Bmat*matrix(ydata,dim(Bmat)[1],dim(Bmat)[2]),2,sum)+xdata[lags+it,] %*% t(Bx) + shocks[it, ]
  }
  if(!is.null(dimnames(y0))){
    dimnames(yhat) <- list(NULL,dimnames(y0)[[2]])
  }
  if(is.ts(y0)){
    yhat <- ts(yhat,start=tsp(y0)[1],freq=tsp(y0)[3])
  }
  return(yhat)
}
