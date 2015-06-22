################################
#Load libraries
################################


################################
#Function
################################
plutusTvv<-function(listData, nit, CBrk =FALSE, seedx=NULL,  seedAcountry= NULL, seedA=NULL, seedH=NULL,H1=10^-2, H2=10^-4, lcA0=NULL, lcLmd=NULL, nSteps=60, nLags=5,verbose=TRUE, nonorm=TRUE, const=FALSE){
	##Function: finds the A0, x, lmd that maximizes the posterior. Allows for both time, and country varying volatility.
		##ListData: ouput from Tvvdata, which stacks country data, and creates a Bx matrix of country dummies
		##nit: nit for csminwelNew
		##CBrk: set to TRUE if you have country breaks
		##seedx: seed for x (concatenated vector of A and lmd values)
		##seedAcountry: country whose residuals will be used to compute the seed A. Must input country abbreviation as string.  
		##H1: Sets scale of A parameter diagonals in inverse Hessian
		##H2: Sets scale of lmd parameter diagonals in inverse Hessian
		##lcA0: logical matrix for A0. If NULL, the default is lower triangular.
		##lcLmd: logical matrix for lmd matrix. If NULL, all TRUE (variances vary in every period)
		##nSteps: number of steps over which to calculate the impulse response output 
		##nLags: number of lags
		##verbose: if TRUE, get the full output from the optimization
		##nonorm: don't normalize with prior	
	##note regarding output: 
		##This program sets Tsigbrk to start at 0.
		##This program does *not* set breaks to start at 0. rfvar3 does this.
		
	
	
################################
#Define
################################
##
#Dates
##
	Tsigbrk=listData$Tsigbrk
	startdate=1 #Tvvdata forces time series to start at 1
	enddate = dim(listData$Y)[1]

##
#Priors
##	
	nVar = dim(listData$Y)[2]
	varnames = dimnames(listData$Y)[[2]]
	
	mnprior = list(tight = 3, decay = .5)
	vprior = list(sig = rep(.01,nVar), w = 0)
	
	asig <- 2
	asd <- outer(vprior$sig, 1/vprior$sig)
		
	urprior <- list(lambda=5, mu=1)
	sigfix <- diag(vprior$sig^2)
	
	names(vprior$sig) = varnames
	dimnames(sigfix) = list(varnames, varnames)
	
	pparams=list(urprior = urprior, asig = asig, mnprior = mnprior, vprior = vprior, asd = asd)



################################
#Sow Seeds
################################

	vars=dimnames(listData$Y)[[2]]
	nVar = dim(listData$Y)[2]
	
	#A0 has a lower triangular default structure 
	if (is.null(lcA0)) {
		lcA0 = matrix(FALSE,nVar,nVar)
		lcA0[lower.tri(lcA0)] = TRUE #lower triangular default structure
	}

###
#Seed A, lmd matrix
##			

#Our operating assumption is that heteroskedastity should idenitfy our A0 matrix (i.e. our model is good; heteroskedasticity is due to contemporaneous relationships between the variables).
		#A * y_t = A+*Y_t-1 + e_t
		#y_t = invA A+ Y_t-1 + invA e_t
		#Hence, residuals from a reduced form VAR should look like invA*e_t
	
		
#I run rfvar3 on the By and Bx matrix and estimate the tOmega matrix by taking the cross product of residuals (divided by # of observations). My A0 estimate is the cholesky lower triangular matrix of tOmega inverse with the diagnoals normalized to one.


#Set Tsigbrk = c(0, Tsigbrk) 
	Tsigbrk0=listData$Tsigbrk
	Tsigbrk=c(0, Tsigbrk0)
	
#Set Countrybrk=Countrybrk 
	if (CBrk==TRUE){
	Countrybrk=listData$Countrybrk
	}
	
#Set Regimes	
	#n=nCountries+1
	T = dim(listData$Y)[1]
	regimes <- c(Tsigbrk, T)
	regimes <- diff(regimes) - nLags
	#n=nCountries+2
	regimes <- c(0, cumsum(regimes))  
		#rfvar3 adds six dummy observations to the end of the data matrix 
	#n=nCountries+1
	nRegimes = length(regimes) - 1
	seedLmd = matrix(0, nVar, nRegimes)

#Create seed A inv
			
	#rfvar3 sets breaks=c(0, Tsigbrk, T)		
	if (CBrk==TRUE){
		seedrfmodel = rfvar3(ydata = listData$Y, xdata=listData$X, lags = nLags, lambda = pparams$urprior$lambda, mu = pparams$urprior$mu, const=const, breaks=Countrybrk)
	}
	else{
		seedrfmodel = rfvar3(ydata = listData$Y, xdata=listData$X, lags = nLags, lambda = pparams$urprior$lambda, mu = pparams$urprior$mu, const=const, breaks=Tsigbrk)
	}
	u=seedrfmodel$u
	seedu=u
	
	if(is.null(seedx)){
		
	if (is.null(seedA)){	
		if (!is.null(seedAcountry)){			
					k= which(seedAcountry==listData$countries, arr.ind=TRUE)
					seedu= seedrfmodel$u[(regimes[k]+1): (regimes[k+1]),]	
					}
								
	#construct temporary Omega							
		tOmega= crossprod(seedu) / dim(seedu)[1]
	#construct temporary A which is the lower cholesky decomposition of the inverse of tOmega	
		tA= t(chol(solve(tOmega)))
	#take the inverse of A. Since the inverse of lower triangular is lower triangular, tAinv will be lower triangular
		tAinv=solve(tA)	
	#take the diagonals of tAinv, as our standard errors
		fullLambda=diag(tAinv)
	#get the lmds by calculating -log(Lambda)	
		fullLmd=-log(fullLambda)
	#normalize the tAinv by setting the diagonals to ones
		tAinv=tAinv/diag(tAinv)	
	#let our seedA be the inverse of this normalized tAinv	
		seedA=solve(tAinv)
	#enforce restrictions	
		seedA[!lcA0]=0 
	#check the lower triangular structure
		diag(seedA)=1}
	#define seedAinv	
		seedAinv=solve(seedA)
	
	
#Create seed Lmd matrix

	if (nRegimes==1){
		seedLmd=fullLmd
	}
	else {
		for(iRegime in (1:nRegimes)){
			#set start of regime, and end of regime
			#recall: regimes denoted as beginning at end of last regime
				iRange=(regimes[iRegime]+1): (regimes[iRegime+1])
				iU=u[iRange,]
			#reduced form variances
				iOmega=crossprod(iU)/dim(iU)[1] 
			
			#GLS it to scale/get str variances. Recall that the lmd matrix is the -log(variance) matrix. 
			iLambda = seedAinv%*%iOmega%*%t(seedAinv) 
			seedLmd[,iRegime]= -log(diag(iLambda))	
		}
		}
	
	if (is.null(lcLmd)){
		#nRegimes=nCountries
	lcLmd=matrix(TRUE, nVar, nRegimes)
	}


##
#Construct seed x
##
	
	seedx=c(c(seedA[lcA0]), c(seedLmd[lcLmd]))
	
	}

##
#Construct seed H. H is the *inverse* Hessian.
##

	nAparams = sum(lcA0)
	nLmdparams = sum(lcLmd)
	
	if (is.null(seedH)){
	seedH = diag(c(rep(H1, nAparams), rep(H2,nLmdparams)))
	}


################################
# Opt
################################
	if (CBrk ==TRUE){
		opt=csminwelNew(bvarwrap5, seedx, seedH, nit=nit, listData=listData, lcA0=lcA0, lcLmd=lcLmd, Tsigbrk=Tsigbrk, breaks=Countrybrk, pparams=pparams, nLags=nLags, nonorm=nonorm)
	}
	else{
		opt=csminwelNew(bvarwrap5, seedx, seedH, nit=nit, listData=listData, lcA0=lcA0, lcLmd=lcLmd, Tsigbrk=Tsigbrk, breaks=Tsigbrk0, pparams=pparams, nLags=nLags, nonorm=nonorm)
	}
	lh = opt$fh
	x = opt$xh
	H= opt$H

##
#Run bvarwrap5
##

	if (CBrk ==TRUE){
		mlmodel = bvarwrap5(x, verbose = TRUE, listData = listData, nLags = nLags, lcA0 = lcA0, lcLmd = lcLmd, Tsigbrk = Tsigbrk, breaks=Countrybrk, pparams = pparams, nonorm=nonorm, const=const)
	}
	else{
		mlmodel = bvarwrap5(x, verbose = TRUE, listData = listData, nLags = nLags, lcA0 = lcA0, lcLmd = lcLmd, Tsigbrk = Tsigbrk, breaks=Tsigbrk0, pparams = pparams, nonorm=nonorm, const=const)
	}
	A = mlmodel$A
	lmd = mlmodel$llmd #ACTUAL lmd-logs
	lambda = mlmodel$lambda #exponentiated 
	relLambda = lambda / lambda[,1] #variances relative to the first period
	vout = mlmodel$vout
	


##
#Calculate bvarwrap5 
##
	
	nLambdas = dim(lambda)[2]
	ir = array(0, dim = c(nVar, nVar, nSteps, nLambdas))
	
	smatlst=list()
	#scales the impulse responses by shock standard error; we have the lambdas=exponentiated log variances
	for (iLambda in (1:nLambdas)){
		#form of the shocks
		
		#structural matrix: In order to make the residuals E we have to multiply A by the lmds to get unit variances
		smat = solve(A) %*% diag(sqrt(lambda[,iLambda]))
		smatlst[[iLambda]]=smat
		By = vout$var$By
			#nvar x nvar x lags matrix of coefficients on lagged y's
		nLags = dim(By)[3]
		for (iLag in (1:nLags)) {
			By[,,iLag] = solve(A) %*% By[,,iLag]
		}
		tempvar = vout$var
		tempvar$By = By
		ir[,,,iLambda] = impulsdtrf(tempvar, smat = smat, nstep = nSteps)
		dimnames(ir[,,,iLambda])[[1]] = vars
	}
	



################################
#Set Return
################################

if (!verbose){
	return(list(A = A, H=H, x=x, lmd = lmd, lambda = lambda, relLambda = relLambda, vout = vout, ir = ir, Tsigbrk=Tsigbrk, startdate= startdate, enddate= enddate))
} else {
	#return full optimization output
	return(list(opt = opt, seedx=seedx, seedA=seedA, seedLmd=seedLmd, A = A, H=H, x=x, lmd = lmd, lambda = lambda, relLambda = relLambda, vout = vout, ir = ir, Tsigbrk=Tsigbrk, regimes=regimes, startdate = startdate, enddate = enddate, smatlst=smatlst))
}
}