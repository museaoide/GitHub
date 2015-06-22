
################################
#Function I
################################

McmcIr<- function(x, results, nSteps=nSteps){
	#function: impulse response 
	#x=x (proposed x or posterior)
	#results: results from plutusTvv
	#nSteps: number of steps in the IR. For us, the number corresponds to number of quarters
	#...: in our case, other inputs for bvarwrap5()	

listData=results$opt$listData
nLags = results$opt$nLags
lcA0 = results$opt$lcA0
lcLmd = results$opt$lcLmd
Tsigbrk = results$opt$Tsigbrk
pparams = results$opt$pparams

mlmodel= bvarwrap5(x, verbose=TRUE, listData=listData, nLags=nLags, lcA0=lcA0, lcLmd=lcLmd, Tsigbrk=Tsigbrk, pparams=pparams, nonorm=TRUE)

A=mlmodel$A
nVar=dim(A)[1]
lmd = mlmodel$llmd 
lambda = mlmodel$lambda 
relLambda = lambda / lambda[,1]
vout = mlmodel$vout

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

McmcMIr<- function(results, nit, Draws=NULL, Sigma=NULL, nSteps=60){
#function: calculates IRs for mcmc draws (tailored to use with plutusTvv and bvarwrap5() output)
	#results: output from optimization through plutusTvv	
	#nit: number of iterations
	#Draws: mcmc draws. If NULL, function will do draws for you.
	#Sigma: Sigma is variance of mean zero gaussian jumps, which is our draw jump function
	#nSteps:IR horizon
	#verbose: TRUE


if (!is.null(Draws)){
nit=dim(Draws$mcout)[1]}

if (is.null(Draws)){
		Draws=Mtrop3(nit=nit, results=results)
		}
		
	verbose=TRUE
	nVar=dim(results$A)[1]
	lambda=results$lambda
	nLambdas= dim(lambda)[2]
	
	irlst=array(0, dim=c(nVar, nVar, nSteps, nLambdas, nit))
		#resp var, shockvar, nSteps, country, nit
	dimnames(irlst)= list(rownames(results$A), colnames(results$A), 1:nSteps, 1:nLambdas, 1:nit)
	
	
	for (it in 1:nit){
		irlst[,,,,it]=McmcIr(x=Draws[[1]][it,], results, nSteps=nSteps)
	}
	
return(irlst)	
}

################################
#Execute
################################