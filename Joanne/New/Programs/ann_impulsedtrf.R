

#By: nvar x nvar x lags matrix of coefficients on lagged y's.  1st dimension is "equation number"
#Bx:      nvar x nx matrix of coefficients on x's
#Fills first the last dimension, by column
#If you resize it will populate by column

B <- vout$By
neq<-dim(B)[1]
nvar<-dim(B)[2]
lags<-dim(B)[3]
dimnB<-dimnames(B)

if(is.null(smat)){
	if (is.null(order)){
		order<-1:neq
	}
	
	#chol of variance approximation with ordering as given
	smat <- t((pchol)(crossprod(vout$u))/dim(vout$u)[1], order))
}
nshock <-dim(smat)[2]

if (dim(smat)[1] != dim(B)[1]) stop ("B and smat matrix conflict")
response <- array(0, dim=c(neq, nshock, nstep+lags-1))

#neq, nshock, nstep+lags-1, the nLag'th element will be the cholesky decomposition
#enough observations such that you can start.
#you'll delete the first nlags-1 set of zeros

response[,,lags] <-smat 

#transpose an array by permuting its dimension
#interchange the last two dimensions
#now neq, nsteps, nshock
response <-aperm(response, c(1,3,2))

##We're gonna move backwards. The first response is given by the cholesky decomp.
##Now we want to go to the second response.

#rhand side variables number of lags*variabes
irhs <- 1:(lags*nvar)
#left hand side variable stands for the 1:nvar that comes after the first lagged set (lags*nvar)
ilhs <- lags*nvar + (1:nvar)


response <- matrix(response, ncol=nshock)
B <- B[,,seq(from=lags, to=1, by=-1)], reverse time order
B<-matrix(B, nrow=nvar)

#recursively solve for the impulse response

for (it in 1:(nstep-1)){
	
	#matrix multiplication
	response[ilhs,]<-B%*% response[irhs,]
	irhs <- irhs + nvar
	ilhs <- ilhs + nvar
}

dim(response)<- c(nvar, nstep+lags-1, nshock)

#switch again the last two dimensions
#taking *out* the last few lines from the dimensions, insert

response <-aperm(response[,-(1:(lags-1)), , drop=FALSE], c(1,3,2))
dimnames(response)<-list(dimnB[[1]], dimnames(smat)[[2]], NULL)
return(response)

