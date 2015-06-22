fcast <- function(y0,By,Bx,xdata,horiz){
### By: equations x variables x lags
### Bx: equations x nx
### xdata: (lags+horiz) x nx
### y0: lags x nvar.  Must be a matrix, even if univariate.
### produces horiz-step-ahead forecast.
### Value is:
### yhat: (horiz+lags) x nvar
  lags <- dim(y0)[1]
  nvar <- dim(y0)[2]
  nx <- dim(Bx)[2]
  yhat <- matrix(0,horiz+lags,nvar)
  yhat[1:lags,] <- y0
  Bmat <- aperm(By,c(3,2,1))
  Bmat <- Bmat[seq(lags,1,-1),,]
  dim(Bmat) <- c(lags*nvar,nvar)
  for (it in 1:horiz){
    ydata <- yhat[it:(it+lags-1),,drop=FALSE]
    yhat[lags+it,] <- apply(Bmat*matrix(ydata,dim(Bmat)[1],dim(Bmat)[2]),2,sum)+xdata[lags+it,]%*%t(Bx)
  }
  if(!is.null(dimnames(y0))){
    dimnames(yhat) <- list(NULL,dimnames(y0)[[2]])
  }
  if(is.ts(y0)){
    yhat <- ts(yhat,start=tsp(y0)[1],freq=tsp(y0)[3])
  }
  return(yhat)
}
