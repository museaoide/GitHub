


irQuants<-function(results, irdraws0){

#results
#irdraws= [a,b,c,d,e]
#a=resp variable
#b=shockvariable
#c=nSteps
#d=country
#e=number of trials

countries=results$opt$listData$countries
nCountries=length(countries)

nVar=dim(results$opt$listData$Y)[[2]]
irquants=array(0, dim=c(nVar, nVar, nCountries, 60, 2))

for (i in 1:nCountries){
	irdraws=irdraws0[,,,i,]
	
	shockvars = 1:nVar
	responsevars = 1:nVar
	nShockvars = length(shockvars)
	nRespvars = length(responsevars)
	
for (iR in 1:nRespvars) {
	iResptrials = irdraws[iR,,,]
	
	for (iS in 1:nShockvars){
				iShocktrials = iResptrials[iS,,]
				#resulting vector (quantiles (1,2), nSteps)
				irquants[iR, iS, i, , ] = t(apply(iShocktrials,1,quantile,probs = c(.025,.975)))
			}
			
			}
			}	

return(irquants=irquants)
}

