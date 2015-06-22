
################################
#Function
################################
bvarwrap6 <- function(x,  listData, nLags, lcA0 = NULL, lcLmdt, lcLmdc, nSecondaryBrk=nSecondaryBrk, nCountries=nCountries, Tsigbrk = NULL, Countrybrk=Countrybrk, pparams = NULL, nonorm=TRUE, verbose=FALSE, const=FALSE) {

##Function: calculates LH|A,lmd. Modification of older bvarwrap function on file.
##Input:
	##listData: output from tvvData_* program
	##nLags: number of lags in structural VAR
	##lcA0: Logical matrix for A0. If NULL, the default is lower triangular.
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
	A = matrix(0, nVar, nVar)
	varnames = dimnames(listData$Y)[[2]]	
	lmdseq=listData$lmdseq
	lmdendseq=listData$lmdendseq
##
#A0
##
	
	if (is.null(lcA0)){
	lcA0=matrix(TRUE, nVar, nVar)
	diag(lcA0)=FALSE
	}
	
	
	nA = sum(lcA0) #number of non-fixed parameters
	A[lcA0] = x[1:nA]	
	diag(A) = 1
	dimnames(A) = list(varnames, varnames)

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
	
	allh = -.5 * sum(((A / pparams$asd)[lcA0])^2) / (pparams$asig)^2 -.5 * (nVar^2 - (length(A) - nA) * (log(2 * pi) + 2 * log(pparams$asig)))

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
	#First lmdc=rep(1...)
	lmdc[lcLmdc]=x[nA+nLmdt+(1:nLmdc)]
	
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

	vout = SVARhtskdmdd2(ydata=listData$Y, xdata=listData$X, lags = nLags, const = const, lmd=lmd, A0 = A, Tsigbrk = Tsigbrk, Countrybrk= Countrybrk, urprior = pparams$urprior, mnprior = pparams$mnprior,vprior = pparams$vprior, train = 0, nonorm=nonorm)
	
#lh is height of posterior density over A0, lmd, A+ at peak.  w is the height of the marginal posterior for A0, lmd, with A+ integrated out.
	lh = -sum(vout$w)
	attr(lh, 'sigpar') = list(A0 = A, lmd = lmd, Tsigbrk = Tsigbrk, Countrybrk=Countrybrk)
	attr(lh, 'prior') = list(mnprior = pparams$mnprior, vprior = pparams$vprior,
	urprior = pparams$urprior, sigfix = pparams$sigfix)
	attr(lh, 'T') = T
	lplmd = -sum(lmd) - sum(exp(-lmd))
	ev = 0 #can be changed
	
#marginal posterior | lmd, A
	lh = lh + ev - lplmd - allh 
	
	if(verbose) {
        ustd <- vout$var$u
        ulevel <- vout$var$uraw
        return(list(lh=lh, vout=vout, A=A, lambda = exp(-lmd), llmd = lmd, lmdt= lmdt , lmdc= lmdc, u=ulevel,
                    ustd=ustd, asig = pparams$asig))} 
   else {
        return(lh)
    }
	
}


