
##################################FUNCTION: RUN VARS, IRs
#9/6/14

rfvar6 <- function(data, countries, credittype="apa", powerof=1, nSteps=60, smat=TRUE){
	#function: runs reduced form VARs, IRs over a system of five variables (log defl, log rgdp, short term rate, log credit2gdp, log real, property prices), and one exogenous variable (fed funds rate)
	#data: input data source
	#countries: input country abbreviation in "", or vector of country abbreviations
	#credittype: BIS credit series ("nbpa", "aha", "ahu", "ana", "apa", "apu", "bpa", "bpu")
	#power of: set credit2gdp variable to power
	#nSteps: number of steps over which to calculate impulse responses
	#output:
		#vout
		#resplst
		#smat
#create matrix to store data	
	smat <- matrix(0, dim(as.matrix(countries))[1],5)
	colnames(smat)<-c("deflTst$pval", "rgdpTst$pval", "strateTst$pval", "c2gTst$pval", "rrpiTst$pval")	
	rownames(smat)<- countries

#Create list objects to store dates
	Startdates=ls()
	Enddates=ls()		
	
#create list to store list
	voutlst<-list()	
	resplst<-list()

	nCountries=length(countries)
	
for (iCountry in countries){
	

		a<- paste(iCountry, "gdp", sep="_")
		b<- paste(iCountry, "defl", sep="_")
		d<- paste(iCountry, "rate", sep="_")
		f<- paste("rp", iCountry, sep="_")
		
			if(!credittype=="nbpa"){
				c<- paste(iCountry, credittype , sep="_")
				credit2gdp <- log((data[,c]/data[,a])^powerof)}		
				
			if(credittype=="nbpa"){
				credit1=NULL
				credit2=NULL
				
				credit1<-paste("apa", iCountry, sep="_")	
				credit2<-paste("bpa", iCountry, sep="_")
				
				if (!is.null(credit1) & !is.null(credit2)){
				credit2gdp <- log(((data[,g]-data[,h])/data[,a])^powerof)
			}
			}
		
		rgdp <- log(data[,a]/data[,b])
		strate <- data[,d]*100
		defl <- log(data[,b])
		rrpi<-log(data[,f])
		ffr<- (data[, "us_ffr"])
		
		datay<- cbind(defl, rgdp, strate, credit2gdp, rrpi, ffr)		
		datay <- ts(datay, start=c(1955,1), freq=4)
		datay <-trimts(datay)
		
#Start/End date
 		
		start=tsp(datay)[1]
		end=tsp(datay)[2]
		Startdates[[iCountry]]=start
		Enddates[[iCountry]]=end


		datax <- datay[,6]
		datay <- datay[, 1:5]
		
#reduced form VAR		
		vout <-rfvar3(ydata=datay, xdata=datax, lags=5, const=TRUE)
#impulse response
		resp <- impulsdtrf(vout, nstep=nSteps)	
##
#Tsts
##

#DeflTst
	yzrone <- array (1, c(5,5,5))
	yzrone[2:5, 1, ]<-0	

	deflTst <- restrictVAR(vout, yzrone=yzrone)
#RgdpTst
	yzrone <- array(1, c(5,5,5))
	yzrone[1,2, ]<-0
	yzrone[3:5,2, ]<-0
	
	rgdpTst <- restrictVAR(vout, yzrone=yzrone)
#StrateTst
	yzrone <- array(1, c(5,5,5))
	yzrone[1:2,3, ]<-0
	yzrone[4:5, 3, ]<-0
	
	strateTst <- restrictVAR(vout, yzrone=yzrone)
#C2gTst
	yzrone <- array(1, c(5,5,5))
	yzrone[1:3,4, ]<-0
	yzrone[5, 4, ]<-0
	
	c2gTst <- restrictVAR(vout, yzrone=yzrone)
#Rrpi Tst
	yzrone <- array (1, c(5,5,5))
	yzrone[1:4, 5, ]<-0	
	
	rrpiTst <- restrictVAR(vout, yzrone=yzrone)
	

		
		smat [iCountry, ]<- c(deflTst$pval, rgdpTst$pval, strateTst$pval, c2gTst$pval, rrpiTst$pval)
		
		resplst[[iCountry]]<-resp
		voutlst[[iCountry]]<-vout
		
	}

	countrysumm<- list(voutlst, resplst, smat=smat, countries=countries, start=Startdates, end=Enddates)
	return(countrysumm)
	
	}
