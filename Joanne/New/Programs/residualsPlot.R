#################Function: Plot Residuals from reduced form

residualsPlot<- function(rfvar6output, plotvar="rgdp"){
	##Function: This program plots the residuals from a reducedform VAR in which breaks denote when country series end. 
	##Plotvar: The variable whose residuals you want to plot.
	##rfvar6output: Results from rfvar6()

		
countries=rfvar6output[[4]]	
nCountries=length(countries)
arr = c(ceiling(nCountries/6), 6)
par (mfrow=arr, col.lab="black", col.main="black", oma=c(1,3,3,2), mar=c(1, 1, 1, 1), tcl = -.1, mgp=c(0,0,0))


for (iCountry in countries)
{
	resid0<-rfvar6output[[1]][[iCountry]]$u[, plotvar]
	start0=as.numeric(rfvar6output[[5]][[iCountry]])
	end0=as.numeric(rfvar6output[[6]][[iCountry]])
	resid <-ts(resid0, start=start0, end=end0, frequency=4)
	
	title=iCountry
	ts.plot(resid, type="h", main=title, xlab='',ylab='')
	abline(a=0, b=0, h=0)
	}
plotTitle=paste("Reduced Form residuals of", plotvar, sep=" ")
filename=paste("RF_residuals", plotvar, sep="_")
title(plotTitle, outer = TRUE, cex = 1.2)
dev.copy2pdf(file = filename)
}

#################Function: Plot Residuals from reduced form by country
residualsPlotc<- function(rfvar6output, beginfilenamewith="Expanded"){
	##Function: This program plots the residuals from a reducedform VAR in which breaks denote when country series end. 
	##          This program is different from residualsPlot in that it plots the residuals of all variables.
	##Input:
	##rfvar6output: results from rfvar6()

	
countries=rfvar6output[[4]]	
nCountries=length(countries)
nShocks=dim(rfvar6output[[1]][[1]]$u)[2]


for (i in 1:nCountries){
	
	arr = c(nShocks, 1)
	par (mfrow=arr, col.lab="black", col.main="black", oma=c(1,3,3,2), mar=c(1, 1, 1, 1), tcl = -.1, mgp=c(0,0,0))

	iCountry<-countries[i]
	
	for (iShock in 1:nShocks){
		
	resid0<-rfvar6output[[1]][[iCountry]]$u[, iShock]
	start0=as.numeric(rfvar6output[[5]][[iCountry]])
	end0=as.numeric(rfvar6output[[6]][[iCountry]])
	resid <-ts(resid0, start=start0, end=end0, frequency=4)
	
	title=toupper(varnames[iShock])
	ts.plot(resid, type="h", main=title, xlab='',ylab='', col="darkblue")
	abline(a=0, b=0, h=0)
	}
	
filename=paste("RF Residuals of", toupper(iCountry), sep=" ")

plotTitle=paste("RF Residuals of", toupper(iCountry), sep=" ")
filename=paste(beginfilenamewith, "_RF_Residuals_", iCountry, ".pdf" ,sep="")
title(plotTitle, outer = TRUE, cex = 1.2)
dev.copy2pdf(file = filename)
}
}

#################Function: Plot Residuals from SVAR

residualsPlot2<- function(results, beginfilenamewith="Expanded"){
	##Function: This program plots the SVAR residuals. SVAR is allowed to have country breaks but
	##          no time breaks.


u=results$vout$var$u
varnames=dimnames(results$opt$listData$Y)[[2]]
Regimes=results$regimes
End=results$opt$listData$End
Start=results$opt$listData$Start


countries=unique(results$opt$listData$countries)
nCountries=length(countries)

nShocks=dim(results$opt$listData$Y)[2]	

	
for (i in 1:nCountries){
	
	arr = c(nShocks, 1)
	par (mfrow=arr, col.lab="black", col.main="black", oma=c(1,3,3,2), mar=c(1, 1, 1, 1), tcl = -.1, mgp=c(0,0,0))

	iCountry<-countries[i]
	
	for (iShock in 1:nShocks){
	resid0<-u[((Regimes[i]+1): Regimes[i+1]), iShock]
	resid1 <-ts(resid0, end=End[[i]], frequency=4)
	title=toupper(varnames[iShock])
	ts.plot(resid1, type="h", main=title, xlab='',ylab='', col="darkblue")
	abline(a=0, b=0, h=0)
	}
	
	plotTitle=paste("SVAR Residuals of", toupper(iCountry), sep=" ")
	filename=paste(beginfilenamewith, "_SVAR_RESIDUALS_", iCountry, ".pdf" ,sep="")
	title(plotTitle, outer = TRUE, cex = 1.2)
	dev.copy2pdf(file = filename)
	
	}
	
}


#################Function: Plot Residuals from SVAR (part c)

residualsPlot2c<- function(results, beginfilenamewith="Expanded"){

##Function: This program plots the SVAR residuals. SVAR is allowed to have country breaks *and*
##time breaks.
##Input: 
##Results: Output from plutusTvv3
		
u=results$vout$var$u
varnames=dimnames(results$opt$listData$Y)[[2]]
Countrybrk=results$opt$listData$Countrybrk
nLags=5

Regimes=c(0, Countrybrk)
Regimes=diff(Regimes)-nLags
Regimes=c(0, cumsum(Regimes))
End=results$opt$listData$End
Start=results$opt$listData$Start

countries=unique(results$opt$listData$countries)
nCountries=length(countries)

SecondaryBrk=results$opt$listData$SecondaryBrk
nSecondaryBrks=length(SecondaryBrk)+1
lmdseq=results$opt$listData$lmdseq
lmdendseq<-results$opt$listData$lmdendseq

nShocks=dim(results$opt$listData$Y)[2]	

k=0

for (i in 1:nCountries){
	
	arr = c(nShocks, 1)
	par (mfrow=arr, col.lab="black", col.main="black", oma=c(1,3,3,2), mar=c(1, 1, 1, 1), tcl = -.1, mgp=c(0,0,0))

	iCountry<-countries[i]
	
	m=lmdendseq[[i]]-lmdseq[[i]]
	
	for (iShock in 1:nShocks){
	resid0<-u[((Regimes[i]+1): Regimes[i+1]), iShock]
	resid1 <-ts(resid0, start=Start[k+1]+(1+nLags)/4, frequency=4)
	title=toupper(varnames[iShock])
	ts.plot(resid1, type="h", main=title, xlab='',ylab='', col="darkblue")
	abline(a=0, b=0, h=0)
	
	for (j in 1:nSecondaryBrks){
	abline(v=SecondaryBrk[j])	
	}
	}
	
	k=k+m
	
plotTitle=paste("SVAR Residuals of", toupper(iCountry), sep=" ")
filename=paste(beginfilenamewith, "_SVAR_RESIDUALS_", iCountry, ".pdf" ,sep="")
title(plotTitle, outer = TRUE, cex = 1.2)
dev.copy2pdf(file = filename)
	}
}






