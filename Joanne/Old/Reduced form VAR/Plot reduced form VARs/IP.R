


################################
#Function
################################

IP <- function(results, postdraw=TRUE, nTrials=1000, countries, shockvars = NULL, responsevars = NULL, shocksize = 1, nSteps = 40, stddev = NULL, shocknames = NULL, varnames = NULL, scale=TRUE) {
	#Function: plots all responses to all shocks
	#results: output from rfvar6

nCountries=length(countries)
for (i in 1:nCountries){
	iCountry=countries[[i]]

	filename=paste("IR6var", iCountry, sep="_")
	
	ir=results[[2]][[i]]

if (postdraw==TRUE){
	
	filename=paste("IR6var", iCountry, "bounds", sep="_" )
	vout=results[[1]][[i]]
	postdraws<-postdrawDiscard(vout, nTrials)
	nDraws = dim(postdraws$smat)[3]
	irdraws = array(NA, c(dim(ir), nDraws))
	for (iDraw in (1:nDraws)){
		drawReg = list(By = postdraws$By[,,,iDraw], Bx = postdraws$Bx[,iDraw])
		irdraws[,,,iDraw] = impulsdtrf(drawReg, smat = t(postdraws$smat[,,iDraw]),
		nstep = dim(irdraws)[3])
	}
	} 
else {irdraws = NULL}

	# Defaults to all structural shocks, all responses
	nVar = dim(ir)[1]
	if (is.null(shockvars)) {shockvars = 1:nVar}
	if (is.null(responsevars)) {responsevars = 1:nVar}
	
	nShockvars = length(shockvars)
	nRespvars = length(responsevars)
	
	if (is.null(varnames)){varnames = rownames(ir)}

#arrange
	gRange = 1:(nSteps)
	arr = c(nRespvars, nShockvars)
	par(mfrow = arr,col.lab="black",col.main="black", 
	  oma=c(1,3,3,2), mar=c(1,.5,.5,.5), tcl=-0.1, mgp=c(0,0,0))
	  
	axisline = rep(0,(nSteps + 2))

for (iRespvar in responsevars) {
	
		iRespseries = shocksize * ir[iRespvar,,]
		iResptrials = shocksize * irdraws[iRespvar,,,]
		
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
			#Plotting each series
			upper = iRespbounds[,2,iShockvar]
			lower = iRespbounds[,1,iShockvar]
		
			#Plot
			plot(iRespseries[iShockvar, 1:nSteps], ylim = c(yMin, yMax), type = 'l', lwd = 1,
			 xlab = '', ylab = '', yaxt = 'n', xaxt = 'n')
			lines(axisline, type = 'l', lwd = 2)
			
			
			#Add scaled axis
		
			if (scale==TRUE){
			iUstd=results[[1]][[iCountry]]$u
			Omega = crossprod(iUstd) / dim(iUstd)[1]
			stddev = sqrt(Omega[iShockvar, iShockvar])	
			}
			else {stddev =1}
			
			if (!is.null(stddev)) {
					 tickspots = axTicks(2)
					 scaleticks = tickspots / stddev
					 scaleticks = signif(round(scaleticks, 2), 2)
				 }
				 
			axis(side=2, at=tickspots, cex.axis=.75)
				
				
			#Plot bounds	
				if (!is.null(upper)) {
					polygon(c(gRange, rev(gRange)), c(lower[gRange], rev(upper[gRange])),  
					col=rgb(0, 0, 1,0.3))
				}
	
						 
			#Text labels
			responsetitle = paste(iCountry, varnames[iRespvar], sep = '_')
			shocktitle= paste (iCountry, varnames[iShockvar], sep='_')
		
			#Adding variable name and axis on leftmost plot
			if (which(shockvars == iShockvar, arr.ind = TRUE) == 1) {
				mtext(responsetitle, side = 2, line = 2, cex = .5)
			}
			
			#Shock name if appropriate
			if (which(responsevars == iRespvar, arr.ind = TRUE) == 1) {
				mtext(shocktitle, side = 3, line = 0, cex = .5)
			}
					
	}
	
	}

bigtitle = paste('Impulse response of', iCountry, 'over', as.character(nSteps), 'periods', sep = ' ')
bigtitle=paste(bigtitle, ".pdf", sep="")
title(bigtitle, outer = TRUE, cex = 1.2)

dev.copy2pdf(file = filename)
}
}



################################
#Execute
################################
