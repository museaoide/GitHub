#################Function: Examine auto-correlation of residuals

autoCorr<- function(rfvar6output,  plotvar="rgdp"){

	
countries=rfvar6output$countries
nCountries=length(countries)
arr = c(ceiling(nCountries/6), 6)
par (mfrow=arr, col.lab="black", col.main="black", oma=c(1,3,3,2), mar=c(1, 1, 1, 1), tcl = -.1, mgp=c(0,0,0))


for (iCountry in countries){
	resid0<-rfvar6output[[1]][[iCountry]]$u[, plotvar]
	start0=as.numeric(rfvar6output$Start[[iCountry]])
	end0=as.numeric(rfvar6output$End[[iCountry]])
	resid <-ts(resid0, start=1, frequency=4)
	acf0=acf(resid, plot=FALSE)
	acf1=acf0$acf
	title=iCountry
	ts.plot(acf1, type="h", main=title, xlab='',ylab='')
	abline(a=0, b=0, h=0)
}
plotTitle=paste("Autocorrelation in residuals of", plotvar, sep=" ")
filename=paste("RF_AC", plotvar, sep="_")
title(plotTitle, outer = TRUE, cex = 1.2)
dev.copy2pdf(file = filename)
}

#################Function: Examine auto-correlation of residuals. Plot by country.

autoCorr2<- function(results){
	##Function: Plots autocorrelation in the residuals with 90, 95% error bands
	##Inputs:
		##Results: Results from plutusTvv
		##Note: It doesn't adjust the time series for breaks...so the residuals nLags before each secondary break should be shifted forward nLags
		
##
#Define
##		

		u=results$vout$var$u
		#varnames=dimnames(results$opt$listData$Y)[[2]]
		varnames=c("Deflator", "Real GDP", "Credit-to-GDP", "Real Property Prices", "Short Term Rate", "Commodity Prices")
		dimnames(u)[[2]]=varnames
		Countrybrk=results$opt$listData$Countrybrk
		nLags=5
		
		Regimes=c(0, Countrybrk)
		Regimes=diff(Regimes)-nLags
		Regimes=c(0, cumsum(Regimes))
		End=results$opt$listData$End
		Start=results$opt$listData$Start
		
		countries=unique(results$opt$listData$countries)
		nCountries=length(countries)
		
		SecondaryBrk<-results$opt$listData$SecondaryBrk
		lmdseq<-results$opt$listData$lmdseq
		lmdendseq<-results$opt$listData$lmdendseq
		
		nVar=dim(results$opt$listData$Y)[2]	


##
#Plot ACF plots by variable
##

for (iVar in 1:nVar){

varname=varnames[[iVar]]
arr = c(ceiling(nCountries/6), 6)
par (mfrow=arr, col.lab="black", col.main="black", oma=c(1,3,3,2), mar=c(1, 1, 1, 1), tcl = -.1, mgp=c(0,0,0))

k=0
for (i in 1:nCountries){
	iCountry<-countries[i]	
	
	m=lmdendseq[[i]]-lmdseq[[i]]
	
	iStart=Start[k+1]
	resid0<-u[((Regimes[i]+1):Regimes[i+1]), iVar]
	resid <-ts(resid0, start=iStart, frequency=4)
	
	nobs <-length(resid0)
	sigma <-1/sqrt(nobs)
	cintupper95 <- 1.96*sigma
	cintlower95 <- -1.96*sigma
	cintupper90 <- 1.65*sigma
	cintlower90 <- -1.65*sigma
	
	title=iCountry
	acf0=acf(resid, plot=FALSE)
	acf1=acf0$acf
	ts.plot(acf1, type="h", main=title, xlab='',ylab='')
	abline(a=0, b=0, h=0)
	abline(h=cintupper95, type="l", lty=6, col="blue")
	abline(h=cintlower95,type="l", lty=6, col="blue")
	abline(h=cintupper90, type="l", lty=6, col="red")
	abline(h=cintlower90,type="l", lty=6, col="red")
		
	k=k+m
}

plotTitle=paste("Autocorrelation in residuals of", varname, sep=" ")
filename=paste("SVAR_AC", varname, sep="_")
title(plotTitle, outer = TRUE, cex = 1.2)
dev.copy2pdf(file = filename)

}
}

