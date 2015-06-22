library(zoo)

IPTvv <- function(results, errorbands=FALSE, nit=10000, shockvars = NULL, responsevars = NULL, shocksize = 1, nSteps = 40, stddev = NULL, postdraws = NULL, shocknames = NULL, varnames = NULL, scale=TRUE) {
	
##Function: plots all specified responses to all specified shocks by country
##Input:
	##results: results from plutusTvv
	##nit: number of draws (default is 10000)
	##errorbands: logical. True if you want bands.
	##shockvars: default is NULL, which means responses to all shocks will be plotted
	##responsevars: default is NULL, which means all responses to shocks will be plotted
	##shocksize: default=1
	##nSteps: horizon of impulse response
	##stddev=NULL
	##postdraws: default is NULL, which means postdraws will be calculated for you 
	##shocknames: default is NULL
	##varnames: default is NULL
	##scale: default is NULL

countries=results$opt$listData$countries
varnames=dimnames(results$opt$listData$Y)[[2]]
Start=results$opt$listData$Start
End=results$opt$listData$End

ir0=results$ir

if (is.null(postdraws) & errorbands==TRUE){
irdraws0=McmcMIr4(results=results, nit=nit, Sigmascale=.2)
}
else{
	irdraws0=postdraws
}

nCountries=length(countries)
for (i in 1:nCountries){
	iCountry=countries[[i]]
	filename=paste("IRTvv", iCountry,i, sep="_")
	ir=ir0[,,,i]
	
	if (errorbands==TRUE){
	irdraws=irdraws0$irlst[,,,i,]}

	# Defaults to all structural shocks, all responses
	nVar = dim(ir)[1]
	if (is.null(shockvars)) {shockvars = 1:nVar}
	if (is.null(responsevars)) {responsevars = 1:nVar}
	
	nShockvars = length(shockvars)
	nRespvars = length(responsevars)
	
	if (is.null(varnames)){varnames = rownames(irdraws0)}

#arrange
	gRange = 1:(nSteps)
	arr = c(nRespvars, nShockvars)
	par(mfrow = arr,col.lab="black",col.main="black", 
	  oma=c(1,3,3,2), mar=c(1,.5,.5,.5), tcl=-0.1, mgp=c(0,0,0))
	  
	axisline = rep(0,(nSteps + 2))

for (iRespvar in responsevars) {
	
		iRespseries = shocksize * ir[iRespvar,,]
		iResptrials=NULL
		if (errorbands==TRUE){
					iResptrials = shocksize * irdraws[iRespvar,,,]}
		
		# Getting error bands, if appropriate
		if (!is.null(iResptrials) & length(iResptrials > 0)){
			iRespbounds = array(NA, c(nSteps, 2, nShockvars))
			for (iShock in shockvars){
				iShocktrials = iResptrials[iShock,,]
				iShockbounds = t(apply(iShocktrials,1,quantile,probs = c(.025,.975)))
				iRespbounds[,,iShock] = iShockbounds[1:(nSteps),]
			}
		} else {
			iRespbounds = NULL
		}
	
		#Plot scale
		yMax = max(iRespseries, iRespbounds, 0) + .0001 * abs(max(iRespseries, iRespbounds))
		yMin =  min(iRespseries, iRespbounds, 0) - .0001 * abs(min(iRespseries, iRespbounds))
		
		
		
		for (iShockvar in shockvars){
		
			#Plot
			plot(iRespseries[iShockvar, 1:nSteps], ylim = c(yMin, yMax), type = 'l', lwd = 1,
			 xlab = '', ylab = '', yaxt = 'n', xaxt = 'n')
			
			lines(axisline, type = 'l', lwd = 2)
							 
		if (errorbands==TRUE){		
			#Plot bounds	
			upper = iRespbounds[,2,iShockvar]
			lower = iRespbounds[,1,iShockvar]
				if (!is.null(upper)) {
					polygon(c(gRange, rev(gRange)), c(lower[gRange], rev(upper[gRange])),  
					col=rgb(0, 0, 1,0.3))
				}
			 }
			#Text labels
			responsetitle = paste(toupper(varnames[iRespvar]), sep = '_')
			shocktitle= paste (toupper(varnames[iShockvar]), sep='_')
		
			#Adding variable name and axis on leftmost plot
			if (which(shockvars == iShockvar, arr.ind = TRUE) == 1) {
				axis(side=2, cex=.75)
				mtext(responsetitle, side = 2, line = 2, cex = .75)
			}
			
			#Add axis on right most plot, responses scaled to response s.e.
			if (which(shockvars == iShockvar, arr.ind = TRUE) == nShockvars) {
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
return(list(irdraws0= irdraws0))
}

