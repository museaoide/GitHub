
################################
#Function
################################
bvarwrap7 <- function(x,  listData, nLags, lcA0 , lcIndA0, lcIntA0, lcLmdt, lcLmdc, nSecondaryBrk, nCountries, Tsigbrk = NULL, Countrybrk=Countrybrk, pparams = NULL, A0tight, nonorm=TRUE, verbose=FALSE, const=FALSE) {

##Function: calculates LH|A,lmd. Modification of older bvarwrap function on file.
##Input:
	##listData: output from tvvData_* program
	##nLags: number of lags in structural VAR
	##lcIndA0: Logical matrix for Individual countries 
	##lcIntA0: Logical matrix for International countries
	##lcLmdt: Logical matrix for time varying variances. 
	##lcLmdc: Logical matrix for country varying variances. 
	##nSecondaryBrk: number of time breaks (# of regimes-1)
	##nCountries: number of countries
	##Tsigbrk: Tsigbrk
	##pparams: Prior parameters
	##nonorm: if TRUE, use dummy observations but do not normalize posterior to make them a proper prior.
	##verbose: 
##Output: if verbose=TRUE
	##lh
	##vout
	##A
	##lambda (which is exp(-lmd)= exp(-log(variance))=variance)
	##llambda (-log(variance))
	##u: u
	##ustd: standardized residuals
	##pparams$asig: asigma prior 


##
#Def
##
	nVar = dim(listData$Y)[2]
	T = dim(listData$Y)[1]
	varnames = dimnames(listData$Y)[[2]]	
	lmdseq=listData$lmdseq
	lmdendseq=listData$lmdendseq
	nCountries=length(unique(listData$countries))
##
#A0
##	
	nA0=sum(lcA0)
	nIntA=sum(lcIntA0)
	nIndA=sum(lcIndA0)	
	nA=nIntA + nCountries*nIndA
	
	#iAbar is the international mean 	
	iAbar=matrix(0,nVar,nVar)
   diag(iAbar)=1
   iAbar[lcIntA0]=x[1:nIntA]
	

            
	
##
#Priors
##
	
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
	
	
	#Adapted for prior
	allh = -.5 * sum(((iAbar / pparams$asd)[lcIntA0])^2) / (pparams$asig)^2 - sum(log(pparams$asd[lcIntA0])) - .5*(nVar^2 - sum(!lcIntA0)) * (log(2 * pi) + 2 * log(pparams$asig))
	q
	#

##
#Lmdt and Lmdc
##


	
	nLmdt=sum(lcLmdt)
	nLmdc=sum(lcLmdc)
	
	#Since the first columns of Lmdt and Lmdc are set to 1, 
	#dim(lcLmdt) is nVar x nRegimes, 
	#where nRegimes=nSecondaryBrk+1 
	#dim(lcLmdc) is nVar x nCountries

	lcLmdt=cbind(matrix(FALSE, nVar, 1), lcLmdt)
	lcLmdc=cbind(matrix(FALSE, nVar, 1), lcLmdc)
	
	
	lmdt=matrix(1, nVar, nSecondaryBrk+1)
	lmdc=matrix(1, nVar, nCountries)
	
	
	#First lmdt=rep(1...)
	lmdt[lcLmdt]=x[nA+(1:nLmdt)]
	dim(lmdt)=c(nVar, nSecondaryBrk+1)
	#First lmdc=rep(1...)
	lmdc[lcLmdc]=x[nA+nLmdt+(1:nLmdc)]
	dim(lmdc)=c(nVar, nCountries)
	
	lmd=NULL
	
	for (iC in 1:nCountries){
			m=lmdseq[[iC]]
			n=lmdendseq[[iC]]
			lmd=cbind(lmd,lmdt[,(m:(n-1))]+lmdc[,iC])
			}
	  
		
##
#vout
##
#set const=FALSE since our Bx matrix is of country dummies	
	vout = SVARhtskdmdd3(listData=listData, lcIndA0= lcIndA0, lcIntA0= lcIntA0, ydata=listData$Y, xdata=listData$X, lags = nLags, const = const, lmd=lmd, x=x, iAbar = iAbar , Tsigbrk = Tsigbrk, Countrybrk= Countrybrk, urprior = pparams$urprior, mnprior = pparams$mnprior,vprior = pparams$vprior, train = 0, nonorm=nonorm)
	
#lh is height of posterior density over A0, lmd, A+ at peak.  w is the height of the marginal posterior for A0, lmd, with A+ integrated out.
	lh = -sum(vout$w)
	attr(lh, 'sigpar') = list(A0 = A, lmd = lmd, Tsigbrk = Tsigbrk, Countrybrk=Countrybrk)
	attr(lh, 'prior') = list(mnprior = pparams$mnprior, vprior = pparams$vprior,
	urprior = pparams$urprior, sigfix = pparams$sigfix)
	attr(lh, 'T') = T
	
	#log posterior
	lplmd = -sum(lmd) - sum(exp(-lmd))
	
	#A0   
    Int=x[nIntA+ (1:(nCountries*nIndA))]
    lpA0 <- -.5 * sum(Int^2)/A0tight
    
	 ev = 0 #can be changed
	
#marginal posterior | lmd, A
	lh = lh + ev - lplmd - lpA0 - allh 
	
	if(verbose) {
        ustd <- vout$var$u
        ulevel <- vout$var$uraw
        return(list(x=x, lh=lh, vout=vout, A=A, lambda = exp(-lmd), llmd = lmd, lmdt= lmdt , lmdc= lmdc, u=ulevel,
                    ustd=ustd, asig = pparams$asig))} 
   else {
        return(lh)
    }
	
}


