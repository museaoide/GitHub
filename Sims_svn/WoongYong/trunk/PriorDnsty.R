PriorDnsty <- function(param,prior.param){

	prd <- vector("numeric",length=21); names(prd) <- names(param[1:21])
	
	prd["phi0"] 		<- dlogitbeta(x=param["phi0"],alpha=prior.param$phi0["alpha"],beta=prior.param$phi0["beta"],location=prior.param$phi0["location"],scale=prior.param$phi0["scale"],log=TRUE)
	prd["phi1"] 		<- dnorm(x=param["phi1"],mean=prior.param$phi1["mu"],sd=prior.param$phi1["sd"],log=TRUE)
	prd["phi2"] 		<- dnorm(x=param["phi2"],mean=prior.param$phi2["mu"],sd=prior.param$phi2["sd"],log=TRUE)	
	prd["gamm"] 		<- dnorm(x=param["gamm"],mean=prior.param$gamm["mu"],sd=prior.param$gamm["sd"],log=TRUE)
	prd["rho"]  		<- dgamma(x=param["rho"],shape=prior.param$rho["shape"],rate=prior.param$rho["rate"],log=TRUE)
	prd["sig"]  		<- dgamma(x=param["sig"],shape=prior.param$sig["shape"],rate=prior.param$sig["rate"],log=TRUE)                           
	prd["bet"]  		<- dgamma(x=param["bet"],shape=prior.param$bet["shape"],rate=prior.param$bet["rate"],log=TRUE)
	prd["delt"]  		<- dgamma(x=param["delt"],shape=prior.param$delt["shape"],rate=prior.param$delt["rate"],log=TRUE)
	prd["omega"] 		<- dbeta(x=param["omega"],shape1=prior.param$omega["shape1"],shape2=prior.param$omega["shape2"],ncp=0,log=TRUE)
  prd["psi"]  		<- dgamma(x=param["psi"],shape=prior.param$psi["shape"],rate=prior.param$psi["rate"],log=TRUE)                      
	prd["pibar"] 		<- dnorm(x=param["pibar"],mean=prior.param$pibar["mu"],sd=prior.param$pibar["sd"],log=TRUE)	
	prd["cbar"] 		<- dnorm(x=param["cbar"],mean=prior.param$cbar["mu"],sd=prior.param$cbar["sd"],log=TRUE)	
	prd["taubar"] 	<- dnorm(x=param["taubar"],mean=prior.param$taubar["mu"],sd=prior.param$taubar["sd"],log=TRUE)	
  prd["rhor"]  		<- dgamma(x=param["rhor"],shape=prior.param$rhor["shape"],rate=prior.param$rhor["rate"],log=TRUE)  
  prd["rhopc"]  	<- dgamma(x=param["rhopc"],shape=prior.param$rhopc["shape"],rate=prior.param$rhopc["rate"],log=TRUE)  
	prd["sig2m"] 		<- dgamma(x=1/param["sig2m"],shape=prior.param$sig2m["shape"],scale=prior.param$sig2m["rate"],log=TRUE) - 2*log(param["sig2m"])
	prd["sig2tau"]	<- dgamma(x=1/param["sig2tau"],shape=prior.param$sig2tau["shape"],scale=prior.param$sig2tau["rate"],log=TRUE) - 2*log(param["sig2tau"])
	prd["sig2r"]		<- dgamma(x=1/param["sig2r"],shape=prior.param$sig2r["shape"],scale=prior.param$sig2r["rate"],log=TRUE) - 2*log(param["sig2r"])
  prd["sig2pc"]		<- dgamma(x=1/param["sig2pc"],shape=prior.param$sig2pc["shape"],scale=prior.param$sig2pc["rate"],log=TRUE) - 2*log(param["sig2pc"])
	#prd["sig2pc"]		<- do.call(dgamma,list(x=1/param["sig2pc"],shape=prior.param$sig2pc["shape"],rate=prior.param$sig2pc["rate"],log=TRUE)) - 2*log(param["sig2pc"])
	prd["rhob"]			<- dbeta(x=param["rhob"],shape1=prior.param$rhob["shape1"],shape2=prior.param$rhob["shape2"],ncp=0,log=TRUE)
	prd["sig2b"]		<- dgamma(x=1/param["sig2b"],shape=prior.param$sig2b["shape"],scale=prior.param$sig2b["rate"],log=TRUE) - 2*log(param["sig2b"])
	
	return(prd)
}