tracePlot2c<- function(postdraws){
	nVars=dim(postdraws$draws$mcout)[2]
	draws=postdraws$draws$mcout
	for (i in 1:nVars){
		arr = c(2, 1)
		par (mfrow=arr, col.lab="black", col.main="black", oma=c(1,3,3,2), mar=c(1, 1, 1, 1), tcl = -.1, mgp=c(0,0,0))
		title1=paste("Parameter", i, "Traceplot and ACF Plot", sep=" ")
		plot(draws[,i], main=title1)
		acf(draws[,i])
		filename=paste("Parameter_", i, ".pdf", sep="")
		dev.copy2pdf(file = filename)
	}
	
}