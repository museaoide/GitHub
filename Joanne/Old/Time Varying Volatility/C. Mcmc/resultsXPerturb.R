resultsXPerturb<- function(results, Sigmascale=3){

		n=length(results$opt$x)
		results$opt$xh=results$opt$x+ rmvnorm(1,mean=rep(0,n),sigma=Sigmascale*results$opt$H)
		results$opt$xh=results$opt$x[1,]
		
		x0=results$opt$xh
		listData = results$opt$listData
		nLags = results$opt$nLags
		lcA0 = results$opt$lcA0
		lcLmdt = results$opt$lcLmdt
		lcLmdc=results$opt$lcLmdc
		Tsigbrk = results$opt$Tsigbrk
		Countrybrk = results$opt$Countrybrk
		
		SecondaryBrk=results$opt$listData$SecondaryBrk
		nSecondaryBrk=length(SecondaryBrk)
		
		countries=results$opt$listData$countries
		nCountries=length(unique(countries))

		
		pparams = results$opt$pparams
		lhfcn=bvarwrap6

		results$opt$fh=lhfcn(x0, listData=listData, nLags=nLags, lcA0=lcA0, lcLmdt=lcLmdt, lcLmd=lcLmdc, Tsigbrk=Tsigbrk, Countrybrk= Countrybrk , nSecondaryBrk=nSecondaryBrk, nCountries = nCountries ,pparams=pparams, nonorm=TRUE)

return(results)
}

