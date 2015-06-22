#PLOT rGDP 


rawPlot3<- function(data, countries, credittype="apa", powerof=1){
	#Function: plot raw data. Modified to plot general government public sector debt (raw and logged)
	#Data
	#Countries: Vector of country abbreviations

for (i in 1:nCountries)
{
	iCountry=countries[i]
	a<- paste(iCountry, "gdp", sep="_")
	b<- paste(iCountry, "defl", sep="_")
	c<- paste(iCountry, credittype , sep="_")
	d<- paste(iCountry, "rate", sep="_")
	k<- paste(iCountry, "gd2gdp", sep="_")
	e<-paste("rp", iCountry, sep="_")
	
	ngdp <- data[,a]
	rgdp <- data[,a]/data[,b]
	credit <- data[,c]
	credit2gdp<-data[,c]/data[,a]
	gd2gdp<-data[,k]/100
	defl <-data[,b]
	rrpi <- data[,e]
	strate <- data[,d]*100
	logrgdp <- log(data[,a]/data[,b])
	logcredit<-log(data[,c])
	logcredit2gdp <- log((data[,c]/data[,a])^powerof)
	loggd2gdp <- log(gd2gdp)
	logdefl <- log(data[,b])
	logrrpi <- log(data[,e])
	datanew <- cbind(ngdp, rgdp, credit, credit2gdp, gd2gdp, defl, strate, rrpi, logrgdp, logcredit, logcredit2gdp, loggd2gdp, logdefl, logrrpi, strate)
	datanew <- ts(datanew, start= c(1954,4), freq=4)
	datanew<-trimts(datanew)
	
	nVars=8
	nLogVars=7
	varnames=c("Nominal GDP", "Real GDP", "Credit", "Credit-to-Gdp", "General Government Public Sector Debt", "Deflator", "Short-term rate", "Real Property Prices", "Log Real GDP", "Log Credit", "Log Credit-to-Gdp", "Log Government Debt-to-GDP","Log Deflator", "Log Real Estate Prices", "Short-term rate")
	
	
	arr = c(nVars, 1)
	par (mfrow=arr, col.lab="black", col.main="black", oma=c(1,3,3,2), mar=c(1, 1, 1, 1), tcl = -.1, mgp=c(0,0,0))

	for (iVar in 1:nVars){
	title=varnames[iVar]
	ts.plot(datanew[,iVar], type="l", main=title, xlab='',ylab='')
	abline(a=0, b=0, h=0)
	abline(v=c(2006,1))
	abline(v=c(2008,1))
	}

filename=paste("Plot of", toupper(iCountry), "Data", sep=" ")
plotTitle=paste("Plot", toupper(iCountry),  "Data", sep=" ")
filename=paste(iCountry, "_RAW", ".pdf", sep="")
title(plotTitle, outer = TRUE, cex = 1.2)
dev.copy2pdf(file = filename)

	for (iVar in ((nVars+1):(nVars+nLogVars))){
	title= varnames[iVar]
	ts.plot(datanew[,iVar], type="l", main=title, xlab='',ylab='')
	abline(a=0, b=0, h=0)
	abline(v=c(2006,1))
	abline(v=c(2008,1))
	}
filename=paste("Plot of", toupper(iCountry),  "SVAR Data",  sep=" ")
plotTitle=paste("Plot", toupper(iCountry),  "SVAR Data", sep=" ")
filename=paste(iCountry, "_SVARDATA", ".pdf" ,sep="")
title(plotTitle, outer = TRUE, cex = 1.2)
dev.copy2pdf(file = filename)

}
}

##################################EXECUTE

countries6e2_16=c("aus","be","can","dan","fr","it","jpn","nl","no","pt","se","uk","us")
rawPlot3(data, countries6e2_16)