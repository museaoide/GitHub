
library(zoo)

IPTvv2 <- function(results, irquants, nSteps = 60) {
	
countries=results$opt$listData$countries
varnames=dimnames(results$opt$listData$Y)[[2]]
Start=results$opt$listData$Start
End=results$opt$listData$End
ir0=results$ir
nCountries=length(countries)

for (i in 1:nCountries){
	iCountry=countries[[i]]
	filename=paste("IRTvv", iCountry, sep="_")
	ir=ir0[,,,i]

	# Defaults to all structural shocks, all responses
	nVar = dim(ir)[1]
	shockvars = 1:nVar
	responsevars = 1:nVar
	nShockvars = length(shockvars)
	nRespvars = length(responsevars)


#arrange
	gRange = 1:(nSteps)
	arr = c(nRespvars, nShockvars)
	par(mfrow = arr,col.lab="black",col.main="black", 
	  oma=c(1,3,3,2), mar=c(1,.5,.5,.5), tcl=-0.1, mgp=c(0,0,0))
	  
	axisline = rep(0,(nSteps + 2))

for (iR in 1:nRespvars) {
	
		iRespseries = ir[iR,,]
		iResptrials=NULL
				
		# Getting error bands, if appropriate
		
			iRespbounds = array(NA, c(nSteps, 2, nShockvars))
			for (iS in 1: nShockvars){
				#[,,,nSteps, 1:2]
				iSbounds = irquants[iR, iS, i, ,1:2]
				iRespbounds[,,iS] = iSbounds[1:(nSteps),]
			}
		 
		
		
		#Plot scale
		yMax = max(iRespseries, iRespbounds, 0) + .0001 * abs(max(iRespseries, iRespbounds))
		yMin =  min(iRespseries, iRespbounds, 0) - .0001 * abs(min(iRespseries, iRespbounds))
		
		
		for (iS in shockvars){
		
			#Plot
			plot(iRespseries[iS, 1:nSteps], ylim = c(yMin, yMax), type = 'l', lwd = 1,
			 xlab = '', ylab = '', yaxt = 'n', xaxt = 'n')
			
			lines(axisline, type = 'l', lwd = 2)
							 
			#Plot bounds	
			upper = iRespbounds[,2,iS]
			lower = iRespbounds[,1,iS]
				if (!is.null(upper)) {
					polygon(c(gRange, rev(gRange)), c(lower[gRange], rev(upper[gRange])),  
					col=rgb(0, 0, 1,0.3))
				}
			
			#Text labels
			responsetitle = paste(toupper(varnames[iR]), sep = '_')
			shocktitle= paste(toupper(varnames[iS]), sep='_')
		
			#Adding variable name and axis on leftmost plot
			if (which(shockvars == iS, arr.ind = TRUE) == 1) {
				axis(side=2, cex=.75)
				mtext(responsetitle, side = 2, line = 2, cex = .75)
			}
			
			#Add axis on right most plot, responses scaled to response s.e.
			if (which(shockvars == iS, arr.ind = TRUE) == nShockvars) {
				axis(side=4, cex.axis=.75)
			}
					
	}
	
	}

dateRange=paste("(", as.yearqtr(Start[i]), "-", as.yearqtr(End[i]), ")", sep="")

bigtitle = paste("Impulse response of ", toupper(iCountry), dateRange, " over ", as.character(nSteps), " quarters", sep ="")
filename=paste(toupper(iCountry), dateRange, ".pdf", sep='')
title(bigtitle, outer = TRUE, cex = 1.2)

dev.copy2pdf(file = filename)
}
return(list(irquants= irquants))
}

