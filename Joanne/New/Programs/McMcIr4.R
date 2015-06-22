
################################
#Function I
################################

McmcIr4<- function(x, results, nSteps=nSteps){
	##Function: Implements Metropolis Hastings algorithm to obtain posterior draws 
	##          for use in generating posterior impulse response error bands
	##Input:
	##x=x (proposed x or posterior mode)
	##results: results from plutusTvv
	##nSteps: number of steps in the IR. For us, the number corresponds to number of quarters
	##...: in our case, other inputs for bvarwrap6()	
	
##
#Extract values from optimization result
##	
		listData=results$opt$listData
		nLags = results$opt$nLags
		lcA0 = results$opt$lcA0
		lcLmdt = results$opt$lcLmdt
		lcLmdc = results$opt$lcLmdc
		nSecondaryBrk = results$opt$nSecondaryBrk
		nCountries=results$opt$nCountries
		Tsigbrk = results$opt$Tsigbrk
		Countrybrk = results$opt$listData$Countrybrk
		pparams = results$opt$pparams


		mlmodel= bvarwrap6(x, verbose=TRUE, listData=listData, nLags=nLags, lcLmdt= lcLmdt, lcLmdc= lcLmdc, nSecondaryBrk=nSecondaryBrk, nCountries=nCountries, Tsigbrk=Tsigbrk, Countrybrk=Countrybrk, pparams=pparams, nonorm=TRUE)
		
		
		A=mlmodel$A
		nVar=dim(A)[1]
		lmd = mlmodel$llmd 
		lambda = mlmodel$lambda 
		relLambda = lambda / lambda[,1]
		vout = mlmodel$vout
		
##
#Impulse response
##		
		nLambdas= dim(lambda)[2]
		ir = array(0, dim=c(nVar, nVar, nSteps, nLambdas))
		dimnames(ir)=list(rownames(A), colnames(A), 1:nSteps, 1:nLambdas)
		
		for (iLambda in (1:nLambdas)){
		##
		#Structural matrix scaled by shock size = Lambda	
			smat = solve (A) %*% diag(sqrt(lambda[, iLambda]))
			By=vout$var$By
			nLags=dim(By)[3]
			for (iLag in (1:nLags)){
				By[,,iLag]=solve(A) %*% By[,,iLag]
			}
			tempvar=vout$var
			tempvar$By=By
			ir[,,,iLambda]=impulsdtrf(tempvar, smat=smat, nstep=nSteps)
			dimnames(ir[,,,iLambda])[[1]]=colnames(mlmodel$A)
		}
		return(ir)
		}
		

################################
#Function II
################################

McmcMIr4<- function(results, nit=10000, Draws=NULL, Sigmascale=.5, nSteps=60, Sigmatype=0){
	##Function: calculates IRs for mcmc draws (tailored to use with plutusTvv and bvarwrap6() output)
		#results: output from optimization through plutusTvv	
		#nit: number of iterations
		#Draws: mcmc draws. If NULL, function will do draws for you.
		#Sigma: Sigma scale (factor to multiply inverse hessian by to get the variance of our draw jump function)
		#nSteps:IR horizon
		#verbose: TRUE


if (is.null(Draws)){
		Draws=Mtrop4(nit=nit,results=results, Sigmascale=Sigmascale, Sigmatype=Sigmatype)
		Draws=Draws$mcout
		}
		
	enit=dim(Draws)[1]	

	verbose=TRUE
	nVar=dim(results$A)[1]
	lambda=results$lambda
	nLambdas= dim(lambda)[2]
	
	irlst=array(0, dim=c(nVar, nVar, nSteps, nLambdas, enit))
		#resp var, shockvar, nSteps, country, nit
	dimnames(irlst)= list(rownames(results$A), colnames(results$A), 1:nSteps, 1:nLambdas, 1:enit)
	
	
	for (it in 1:enit){
		irlst[,,,,it]=McmcIr4(x=Draws[it,], results, nSteps=nSteps)
	}
	
return(list(irlst=irlst, draws=Draws))	
}

################################
#Execute
################################