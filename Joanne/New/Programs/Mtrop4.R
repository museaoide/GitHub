Mtrop4<- function(results, x=NULL,nit=nit, Sigma=NULL, Sigmascale=.25, Sigmatype=0) {
	
	
##Function: Implements Metropolis-Hasting algorithm to get mcmc draws from the posterior 
##Input:
	##nit: number of iterations/draws
	##Sigma: Covariance matrix for the Guassian jumps. If Sigma is null, we set it equal to the inverse Hessian. Recall that the H in csminwelNew is the *inverse* Hessian
	##Sigmascale: Scale of default sigma (ex: if .5, Sigma= .5*invH)
	##Sigmatype:
	##	0 Sets Sigma to inverse hessian * Sigmascale. Requires Sigma to be NULL.
	##	1 Set Sigma to diagonal matrix whose diagonal entries are the diagonals of the inverse hessian * Sigmascale. Requires Sigma to be NULL.
  	##results: PlutusTvv output
##Output:
  		##mcout: mcmc draws
  		##lhlist: likelihood of each draw



##
#Define
##
	lhfcn=bvarwrap6
	verbose=FALSE
	
##	
#Pull from plutusTvv results
##
	if (is.null(x)){
	x0<-results$opt$xh}
	else{x0<-x}
	listData = results$opt$listData
	nLags = results$opt$nLags
	lcA0 = results$opt$lcA0
	lcLmdt = results$opt$lcLmdt
	lcLmdc = results$opt$lcLmdc
	nSecondaryBrk = results$opt$nSecondaryBrk
	nCountries=results$opt$nCountries
	Tsigbrk = results$opt$Tsigbrk
	Countrybrk =results$opt$listData$Countrybrk
	pparams = results$opt$pparams

##
#Generate jumps centered at 0 with covariance Sigma
##
	nx0=length(x0)		
	mcout <- matrix(0, nit, nx0)
	mcout[1, ] <- x0

if (is.null(Sigma)){  
	if (Sigmatype==1){
	H=results$opt$H
	Sigma=matrix(0, dim(H)[1], dim(H)[2])
	diag(Sigma)=diag(H)*Sigmascale
	}
	if (Sigmatype==0){
	H=results$opt$H
	Sigma=Sigmascale*H
	} 
	}
	 
  jumps <- rmvnorm(nit, rep(0, nx0) , Sigma)
  lhlist <- vector("numeric", nit)
  lh0 <- lhfcn(x0, listData=listData, nLags=nLags, lcA0=lcA0, lcLmdt= lcLmdt, lcLmdc=lcLmdc,nSecondaryBrk=nSecondaryBrk, nCountries=nCountries,  Tsigbrk=Tsigbrk, Countrybrk=Countrybrk, pparams=pparams, nonorm=TRUE)
  lhlist[1] <- lh0
  
  for ( it in 2:nit) {
    ## x1 <- drawjump(x0)
    x1 <-x0 + jumps[it, ]
    lh1 <- lhfcn(x1, listData=listData, nLags=nLags, lcA0=lcA0, lcLmdt= lcLmdt, lcLmdc= lcLmdc,nSecondaryBrk=nSecondaryBrk, nCountries=nCountries, Tsigbrk=Tsigbrk, Countrybrk =Countrybrk, pparams=pparams, nonorm=TRUE) 
    	#runif=random uniform (default= draw in [0,1])
	   if ((lh0-lh1>log(runif(1)))){
	      lhlist[it] <- lh1
	      lh0 <- lh1
	      x0 <- x1
	      mcout[it, ] <- x1
	    } else {
	      lhlist[it] <-  lh0
	      mcout[it, ] <- x0
	    }
	    
  }
  
  return(list(mcout=mcout, lhlist=lhlist))
}
