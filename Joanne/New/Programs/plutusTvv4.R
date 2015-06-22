 ################################
#Load libraries
################################
#load(abind)

################################
#Function
################################
plutusTvv4<- function(listData, nit, seedx=NULL,  seedIntA0=NULL, seedH=NULL, seedlmdt=NULL, seedlmdc=NULL, H1=10^-2, H2=10^-4, lcIndA0=NULL, lcIntA0=NULL, lcLmdt=NULL, lcLmdc=NULL, nSteps=60, nLags=5,verbose=TRUE, nonorm=TRUE, const=FALSE, mnprior = list(tight = 3, decay = .5), vprior = list(sig = rep(.01,nVar), w = 0), urprior=list(lambda=5, mu=1), A0tight=.01){
	
##Function: finds the A0, x, lmd that maximizes the posterior using Newtonian optimization (csminwelNew). 
##This function allows for both time, and country varying volatility. 
##This program differs from plutusTvv in that Lmd= Lmdc + Lmdt. (Adding prevents weird discontinuities in the objective function.) 	
	##ListData: ouput from the tvvData_* programs, which prepares data for use in plutusTvv_* programs
	##nit: maximum nit for csminwelNew
	##seedx: seed for x (concatenated vector of A and lmd values)
	##H1: Sets scale of A parameter diagonals in inverse Hessian
	##H2: Sets scale of lmd parameter diagonals in inverse Hessian
	##lcIntA0: Logical matrix for the international A0. Elements that are free to vary. Default is set to TRUE for all off diagonal elements.
	##lcA0: Logical matrix for A0. Elements that are free to vary. Default is TRUE for all off diagonal elements.
	##lcLmdt: Logical matrix for time varying variances. First time regime is fixed to be one. Dimensions are nVars x nRegimes
	##lcLmdc: Logical matrix for country varying variances. First country is fixed to be one. Dimension are nVars x nCountries
	##nSteps: number of steps over which to calculate the impulse response output. Default is 60.
	##nLags: number of lags for the endogenous variables.
	##verbose: if TRUE, get the full output from the optimization
	##nonorm: don't normalize with prior	
	##mnprior
	##vprior
	##urprior: set NULL
	##A0tight:  tightness of A0 prior. Prior on the variance of iIndA0 elements which we believe have independent guassian
    ##                distribution with mean 0, and variance A0tight

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

	lcA0=lcIntA0|lcIndA0

###
#Seed A, lmdc, lmdt matrices
##			

	#Set Tsigbrk = c(0, Tsigbrk)
	#Set Countrybrk=Countrybrk (doesn't begin with 0; rfvar3 will set Countrybrk=c(0, Countrybrk,T))
	Tsigbrk=listData$Tsigbrk
	Countrybrk<-listData$Countrybrk
	T = dim(listData$Y)[1]
	Tsigbrk <- c(0, Tsigbrk)


	#Define regimes denotes starting date of a new regime
	regimes <- c(Tsigbrk, T)
	regimes <- diff(regimes) - nLags
	regimes <- c(0, cumsum(regimes))  
	nRegimes = length(regimes) - 1


	#N. parameters
	nA0=sum(lcA0)
	nIntA=sum(lcIntA0)
	nIndA=sum(lcIndA0)
	

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

	#seedx
	if(is.null(seedx)){
	if (is.null(seedIntA0)){	


		seedrfmodel = rfvar3(ydata = listData$Y, xdata=listData$X, lags = nLags, lambda = pparams$urprior$lambda, mu = pparams$urprior$mu, const=const, breaks=Countrybrk)
		u=seedrfmodel$u
		seedu=u				
		tOmega= crossprod(seedu) / dim(seedu)[1]
		tA= t(chol(solve(tOmega)))
		tAinv=solve(tA)	
		tAinv=tAinv/diag(tAinv)	
		seedIntA0=solve(tAinv)
		diag(seedIntA0)=1
		}
		
	#previously when I restricted model to case where lcIntA0 and lcIndA0 are non overlapping
	#seedA=c(seedA[lcIntA0], rep(seedA[lcIndA0], nCountries))	
	
	seedA=c(seedIntA0[lcIntA0], rep(0, nIntA*nCountries))

##
#Construct seed x, which should contain lmds for the dummy
##
	
	seedx=c(seedA,seedLmdc[lcLmdc],seedLmdt[lcLmdt])
	
	}



##
#Construct seed H. H is the *inverse* Hessian. Set it for when lmd has been constructed 
##

	
#Technically, the number of our lmd parameters is sum of lmdc and lmdt	
	nLmdparams =sum(lcLmdt) + sum(lcLmdc)

	if (is.null(seedH)){
	seedH = diag(c(rep(H1, nIntA+nIndA*nCountries), rep(H2,nLmdparams)))
	}

################################
# Opt
################################
	opt=csminwelNew(bvarwrap7, seedx, seedH, nit=nit, listData=listData, lcA0=lcA0, lcIndA0=lcIndA0, lcIntA0= lcIntA0, lcLmdt=lcLmdt, lcLmdc=lcLmdc,  nSecondaryBrk=nSecondaryBrk, nCountries=nCountries, Tsigbrk=Tsigbrk, Countrybrk = Countrybrk, pparams=pparams, A0tight=A0tight, nLags=nLags, nonorm=nonorm, const=const)
	lh = opt$fh
	x <- opt$xh
	H= opt$H
		
##
#Run bvarwrap6
##
	mlmodel = bvarwrap7(x, verbose = TRUE, listData = listData, lcA0 = lcA0, lcIndA0= lcIndA0, lcIntA0= lcIntA0, lcLmdt=lcLmdt, lcLmdc=lcLmdc, nSecondaryBrk=nSecondaryBrk, nCountries=nCountries, Tsigbrk = Tsigbrk, Countrybrk = Countrybrk, pparams = pparams, A0tight= A0tight, nLags=nLags, nonorm=nonorm, const=const)
	
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
	lmdendseq=listData$lmdendseq
	lmdseq=listData$lmdseq
	kLambdas=0	
	
	for (i in 1:nCountries)	{
		
		iA=matrix(1, nVar, nVar)
		iA[lcIntA0]=x[1:nIntA] + x[nIntA+(i-1)*nIndA+1:nIndA]
		mLambdas=lmdendseq[[i]]-lmdseq[[i]]

	for (iLambda in 1:mLambdas)      {
		smat = solve(iA) %*% diag(sqrt(lambda[,kLambdas+iLambda]))
		smatlst[[kLambdas+iLambda]]=smat
		
		
		By = vout$var$By
			#nvar x nvar x lags matrix of coefficients on lagged y's
		nLags = dim(By)[3]
		
		
		for (iLag in (1:nLags)) {
			By[,,iLag] = solve(iA) %*% By[,,iLag]
		}
		
		
		tempvar = vout$var
		tempvar$By = By
		ir[,,,kLambdas+iLambda] = impulsdtrf(tempvar, smat = smat, nstep = nSteps)
		dimnames(ir[,,,kLambdas+iLambda])[[1]] = varnames
	}
	kLambdas=kLambdas+mLambdas
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
