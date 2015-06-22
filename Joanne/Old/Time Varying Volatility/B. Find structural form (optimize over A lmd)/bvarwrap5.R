
################################
#Function
################################
bvarwrap5 <- function(x,  listData, nLags, lcA0 = NULL, lcLmd = NULL,  Tsigbrk = NULL, breaks=breaks, pparams = NULL, nonorm=TRUE, verbose=FALSE, const=FALSE) {

#function: calculates LH|A,lmd. Modification of older bvarwrap function on file.
#listData: output from TvvVar
#nLags: number of lags in structural VAR
#lcA0: logical matrix for A0 matrix
#lcLmd: logical matrix for A0 matrix
#Tsigbrk: Structural breaks
#pparams: Prior parameters
#nonorm: if TRUE, use dummy observations but do not normalize posterior to make them a proper prior.
#verbose: 
#output: if verbose=TRUE
	#lh
	#vout
	#A
	#lambda 
	#llambda (-log(variance))
	#u: u
	#ustd: standardized residuals
	#pparams$asig: asigma prior 


##
#Def
##
	nVar = dim(listData$Y)[2]
	T = dim(listData$Y)[1]
	A = matrix(0, nVar, nVar)
	varnames = dimnames(listData$Y)[[2]]	
	
##
#A0
##
	diag(A) = 1
	if (is.null(lcA0)){
	lcA0=matrix(TRUE, nVar, nVar)}
	nA = sum(lcA0) #number of non-fixed parameters
	A[lcA0] = x[1:nA]	
	dimnames(A) = list(varnames, varnames)

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

allh = -.5 * sum(((A / pparams$asd)[lcA0])^2) / (pparams$asig)^2 -.5 * (nVar^2 - (length(A) - nA) * (log(2 * pi) + 2 * log(pparams$asig)))

##
#lmd
##
	nLmd = sum(lcLmd)
	nSig = length(Tsigbrk)
	lmd = matrix(0, nVar, nSig)
	lmd[lcLmd] = x[nA + (1 : nLmd)]
	
	if (!all(lcLmd)){
		rest = which(!lcLmd, arr.ind = TRUE)
		nRest = dim(rest)[1]
		for (iRest in (1:nRest)) {
			lmd[rest[iRest,1], rest[iRest, 2]] = lmd[rest[iRest,1], rest[iRest, 2] - 1]
		}
	}
	
		
##
#vout
##

#set const=FALSE since our Bx matrix is of country dummies	
	vout = SVARhtskdmdd(ydata=listData$Y, xdata=listData$X, lags = nLags, const = const, A0 = A,
	lmd = lmd, Tsigbrk = Tsigbrk, breaks=breaks, urprior = pparams$urprior, mnprior = pparams$mnprior,
	vprior = pparams$vprior, train = 0, nonorm=nonorm)
		
#lh is height of posterior density over A0, lmd, A+ at peak.  w is the height of the marginal posterior for A0, lmd, with A+ integrated out.
	lh = -sum(vout$w)
	attr(lh, 'sigpar') = list(A0 = A, lmd = lmd, Tsigbrk = Tsigbrk, breaks=breaks)

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
        return(list(lh=lh, vout=vout, A=A, lambda = exp(-lmd), llmd = lmd, u=ulevel,
                    ustd=ustd, asig = pparams$asig))} 
   else {
        return(lh)
    }
	
}


