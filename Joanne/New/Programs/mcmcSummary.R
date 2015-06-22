mcmcSummary<- function(draws, thin=1, burnin=0){
	##Function: outputs the "effective size" of mcmc draws. Thins and deletes draws ("burn-in") if desired.
	##Input:
		##Thin: Thinning factor. 
		##Burnin: Burnin


nit=dim(draws$mcout)[1]
enit=(nit-burnin)/thin

draws$mcout=draws$mcout[-(1:burnin), ]    
draws$mcout=draws$mcout[thin*(1:enit), ]    
draws$lhlist=draws$lhlist[thin*(1:enit)]         
draws$lhlist=draws$lhlist[-(1:burnin)]    

mcmcdraws=draws$mcout
mcmcdraws=as.mcmc(mcmcdraws)
return(list(effectiveSize(mcmcdraws), draws))
}

