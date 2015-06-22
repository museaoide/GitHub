
################################
#Load libraries
################################
library(mvtnorm)

################################
#Function
################################
Mtrop3 <- function( nit=nit, Sigma=NULL, Sigmascale=.5, results) {
	#function: Mcmc draws from posterior 
	#nit: number of iterations/draws
	#Sigma: Covariance matrix for the Guassian jumps. If Sigma is null, we set it equal to the inverse Hessian.Recall that the H in csminwelNew is the *inverse* Hessian
	#Sigmascale: Scale of default sigma (ex: if .5, Sigma= .5*invH)
  	#results: PlutusTvv output
  	#output: list
  		#mcout: mcmc draws
  		#lhlist: likelihood of each draw

#Define
	lhfcn=bvarwrap5
	verbose=FALSE
	
	
#Pull from plutusTvv results
	
	x0=results$opt$xh
	listData = results$opt$listData
	nLags = results$opt$nLags
	lcA0 = results$opt$lcA0
	lcLmd = results$opt$lcLmd
	Tsigbrk = results$opt$Tsigbrk
	pparams = results$opt$pparams


#Generate jumps centered at 0 with covariance Sigma

	nx0=length(x0)		
	mcout <- matrix(0, nit, nx0)
	mcout[1, ] <- x0

if (is.null(Sigma)){  
	Sigma=Sigmascale*results$opt$H
	  } 
  			
  jumps <- rmvnorm(nit, rep(0, nx0) , Sigma)
  lhlist <- vector("numeric", nit)
  lh0 <- lhfcn(x0, listData=listData, nLags=nLags, lcA0=lcA0, lcLmd=lcLmd, Tsigbrk=Tsigbrk, pparams=pparams, nonorm=TRUE)
  lhlist[1] <- lh0
  
  for ( it in 2:nit) {
    ## x1 <- drawjump(x0)
    x1 <-x0 + jumps[it, ]
    lh1 <- lhfcn(x1, listData=listData, nLags=nLags, lcA0=lcA0, lcLmd=lcLmd, Tsigbrk=Tsigbrk, pparams=pparams, nonorm=TRUE)    
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

################################
#Execute
################################

