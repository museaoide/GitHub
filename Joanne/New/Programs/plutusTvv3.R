 ################################
#Load libraries
################################
#load(abind)

################################
#Function
################################
plutusTvv3<- function(listData, nit, seedx=NULL, seedAcountry= NULL, seedA=NULL, seedH=NULL,H1=10^-2, H2=10^-4, lcA0=NULL, lcLmdt=NULL, lcLmdc=NULL, nSteps=60, nLags=5,verbose=TRUE, nonorm=TRUE, const=FALSE, mnprior = list(tight = 3, decay = .5), vprior = list(sig = rep(.01,nVar), w = 0), urprior=list(lambda=5, mu=1)){
	
##Function: finds the A0, x, lmd that maximizes the posterior using Newtonian optimization (csminwelNew). 
##This function allows for both time, and country varying volatility. 
##This program differs from plutusTvv in that Lmd= Lmdc + Lmdt. (Adding prevents weird discontinuities in the objective function.) 	
	##ListData: ouput from the tvvData_* programs, which prepares data for use in plutusTvv_* programs
	##nit: maximum nit for csminwelNew
	##seedx: seed for x (concatenated vector of A and lmd values)
	##seedAcountry: country whose residuals will be used to compute the seed A. Must input country abbreviation as string.  
	##H1: Sets scale of A parameter diagonals in inverse Hessian
	##H2: Sets scale of lmd parameter diagonals in inverse Hessian
	##lcA0: Logical matrix for A0. If NULL, the default is lower triangular.
	##lcLmdt: Logical matrix for time varying variances. First time regime is fixed to be one. Dimensions are nVars x nRegimes
	##lcLmdc: Logical matrix for country varying variances. First country is fixed to be one. Dimension are nVars x nCountries
	##nSteps: number of steps over which to calculate the impulse response output. Default is 60.
	##nLags: number of lags for the endogenous variables.
	##verbose: if TRUE, get the full output from the optimization
	##nonorm: don't normalize with prior	
	##mnprior
	##vprior
	##urprior: set NULL
##Note regarding output: 
	##This program sets *Tsigbrk* to start at 0.
	##This program does *not* set *breaks* to start at 0. rfvar3 does this.
		
	
	
	
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
	
	
	asig <- 2
	asd <- outer(vprior$sig, 1/vprior$sig)
	sigfix <- diag(vprior$sig^2)
	
	names(vprior$sig) = varnames
	dimnames(sigfix) = list(varnames, varnames)
	
	pparams=list(urprior = urprior, asig = asig, mnprior = mnprior, vprior = vprior, asd = asd)



################################
#Sow Seeds
################################

	
	#A0 has free structure
	if (is.null(lcA0)) {
		lcA0 = matrix(TRUE,nVar,nVar)
		diag(lcA0)=FALSE
	}

###
#Seed A, lmdc, lmdt matrices
##			

	#Set Tsigbrk = c(0, Tsigbrk)
	#Set Countrybrk=Countrybrk (doesn't begin with 0; rfvar3 will set Countrybrk=c(0, Countrybrk,T))
	Tsigbrk=listData$Tsigbrk
	Countrybrk=listData$Countrybrk
	T = dim(listData$Y)[1]
	Tsigbrk <- c(0, Tsigbrk)
	Countrybrk <- Countrybrk

	#Define regimes denotes starting date of a new regime
	regimes <- c(Tsigbrk, T)
	regimes <- diff(regimes) - nLags
	regimes <- c(0, cumsum(regimes))  
	
	nRegimes = length(regimes) - 1

	#Create seed A inv	
	seedrfmodel = rfvar3(ydata = listData$Y, xdata=listData$X, lags = nLags, lambda = pparams$urprior$lambda, mu = pparams$urprior$mu, const=const, breaks=Countrybrk)
	u=seedrfmodel$u
	seedu=u

	#Create seed Lmd matrix
	SecondaryBrk=listData$SecondaryBrk
	nSecondaryBrk=length(SecondaryBrk)
		
	countries=listData$countries
	nCountries=length(unique(countries))
	
	#Time lmd. Arbitrary.
		seedLmdt=matrix(rnorm(1, mean=.5, sd=1), nVar, nSecondaryBrk)
		
	#Country lmd. Arbitrary.
		seedLmdc=matrix(rnorm(1.5, mean=2, sd=2), nVar, nCountries-1)
	
	#Set the lcLmdt and lcLmdc logical matrices
	#Since the first lcLmdt is set to one, just set the logicals for the other Lmdts
	#Since the first lcLmdc is set to one, just set the logicals for the other Lmdcs
	if (is.null(lcLmdt)){lcLmdt=matrix(TRUE, nVar, nSecondaryBrk)}
	if (is.null(lcLmdc)){lcLmdc=matrix(TRUE, nVar, nCountries-1)}


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
		
		#set the diagonals to ones
		tAinv=tAinv/diag(tAinv)	
		
		#let our seedA be the inverse of this normalized tAinv	
		seedA=solve(tAinv)
		
		diag(seedA)=1}

##
#Construct seed x, which should contain lmds for the dummy
##
	
	seedx=c(seedA[lcA0],seedLmdc[lcLmdc],seedLmdt[lcLmdt])
	
	}



##
#Construct seed H. H is the *inverse* Hessian. Set it for when lmd has been constructed 
##


	nAparams = sum(lcA0)
	
#Technically, the number of our lmd parameters is sum of lmdc and lmdt	
	nLmdparams =sum(lcLmdt) + sum(lcLmdc)

	if (is.null(seedH)){
	seedH = diag(c(rep(H1, nAparams), rep(H2,nLmdparams)))
	}

################################
# Opt
################################
	opt=csminwelNew(bvarwrap6, seedx, seedH, nit=nit, listData=listData, lcA0=lcA0, lcLmdt=lcLmdt, lcLmdc=lcLmdc,  nSecondaryBrk=nSecondaryBrk, nCountries=nCountries, Tsigbrk=Tsigbrk, Countrybrk = Countrybrk, pparams=pparams, nLags=nLags, nonorm=nonorm)
	lh = opt$fh
	x = opt$xh
	H= opt$H
		
##
#Run bvarwrap6
##
	mlmodel = bvarwrap6(x, verbose = TRUE, listData = listData, nLags = nLags, lcA0 = lcA0, lcLmdt=lcLmdt, lcLmdc=lcLmdc, nSecondaryBrk=nSecondaryBrk, nCountries=nCountries, Tsigbrk = Tsigbrk, Countrybrk = Countrybrk, pparams = pparams, nonorm=nonorm, const=const)
	A = mlmodel$A
	lmd = mlmodel$llmd #ACTUAL lmd-logs
	lambda = mlmodel$lambda #exponentiated 
	relLambda = lambda / lambda[,1] #variances relative to the first period
	vout = mlmodel$vout
	lmdc=mlmodel$lmdc
	lmdt=mlmodel$lmdt
	

	
	nLambdas = dim(lambda)[2]
	ir = array(0, dim = c(nVar, nVar, nSteps, nLambdas))
	
	smatlst=list()
	#scales the impulse responses by shock standard error; we have the lambdas=exponentiated log variances
	for (iLambda in (1:nLambdas)){
		#form of the shocks
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
		dimnames(ir[,,,iLambda])[[1]] = varnames
	}
	



################################
#Set Return
################################

if (!verbose){
	return(list(A = A, H=H, x=x, lmd = lmd, lambda = lambda, relLambda = relLambda, vout = vout, ir = ir, Tsigbrk=Tsigbrk, startdate= startdate, enddate= enddate))
} else {
	#return full optimization output
	return(list(opt = opt, seedx=seedx, seedA=seedA, A = A, H=H, x=x, lmd = lmd, lmdc=lmdc, lmdt=lmdt, lambda = lambda, relLambda = relLambda, vout = vout, ir = ir, Tsigbrk=Tsigbrk, regimes=regimes, startdate = startdate, enddate = enddate, smatlst=smatlst))
}
}


################################
#Execute
################################
