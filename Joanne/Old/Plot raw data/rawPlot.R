#PLOT rGDP 


rawPlot<- function(countries, data, commdata, raw=FALSE, log=TRUE, commIndexWeights =c(1/4, 1/4, 1/4, 1/4),credittype="apa", powerof=1){
	#Function: plot raw data 
	#Data
	#Countries: Vector of country abbreviations
	
nCountries=length(countries)
for (i in 1:nCountries){
	iCountry=countries[i]
	a<- paste(iCountry, "gdp", sep="_")
	b<- paste(iCountry, "defl", sep="_")
	c<- paste(iCountry, credittype , sep="_")
	d<- paste(iCountry, "rate", sep="_")
	e<-paste("rp", iCountry, sep="_")
	
	#exchange rates are USD/LCU
	g<-paste("ex", iCountry, sep="_")
	
	ngdp <- data[,a]
	rgdp <- data[,a]/data[,b]
	credit <- data[,c]
	credit2gdp<-data[,c]/data[,a]
	defl <-data[,b]
	rrpi <- data[,e]
	strate <- data[,d]*100
	logrgdp <- log(data[,a]/data[,b])
	logcredit<-log(data[,c])
	logcredit2gdp <- log((data[,c]/data[,a])^powerof)
	logdefl <- log(data[,b])
	logrrpi <- log(data[,e])

	comm <-commdata[,g]*(commIndexWeights[1]*commdata[,"iENERGY"] + commIndexWeights[2]*commdata[,"iNONFUEL"] + commIndexWeights[3]*commdata[,"iPRECIOUSMET"]+commIndexWeights[4]*commdata[,"Copper"])
	comm <- aggregate(comm, nfreq=4)
#commodityIndex Data (in USD)
	logcomm <- log(commdata[,g]*(commIndexWeights[1]*commdata[,"iENERGY"] + commIndexWeights[2]*commdata[,"iNONFUEL"] + commIndexWeights[3]*commdata[,"iPRECIOUSMET"]+ commIndexWeights[4]*commdata[,"Copper"]))
	logcomm <- aggregate(logcomm, nfreq=4)
		
	datanew <- cbind(rgdp, credit2gdp, defl, strate, rrpi)
	datanew <- ts(datanew, start= c(1954,4), freq=4)
	datanew <- ts.union(datanew, comm)
	datanew<-trimts(datanew)
	
	datanewlog<-cbind(logrgdp, logcredit2gdp, logdefl, strate, logrrpi)
	datanewlog <- ts(datanewlog, start= c(1954,4), freq=4)
	datanewlog<- ts.union(datanewlog, logcomm) 
	datanewlog<-trimts(datanewlog)
	
	
	nVars=dim(datanew)[2]
	nLogVars=dim(datanewlog)[2]
	
	varnames=list("rGDP", "C2G", "Defl", "Rate", "rPP", "CP")
	logvarnames=list("log rGDP", "log C2G", "log Defl", "log Rate", "log rPP", "log CP")
	
if (raw){	
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
}

if (log){
	arr = c(nLogVars, 1)
	par (mfrow=arr, col.lab="black", col.main="black", oma=c(1,3,3,2), mar=c(1, 1, 1, 1), tcl = -.1, mgp=c(0,0,0))

	
	for (iLogVar in 1:nLogVars){
	title= logvarnames[iLogVar]
	ts.plot(datanewlog[, iLogVar], type="l", main=title, xlab='',ylab='')
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
}

##################################EXECUTE

