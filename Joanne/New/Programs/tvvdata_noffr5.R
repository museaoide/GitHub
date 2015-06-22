################################
#Function
################################

tvvData_noffr5<-function(countries, data, commdata, commIndexWeights =c(1/4, 1/4, 1/4, 1/4), freq=4,nLags=5,SecondaryBrk=c(0, 1984, 1998.75, 2007.75, 2015), begin=NULL, terminate=NULL, credittype="apa", powerof=1){
	
##Function: Stacks *differenced* country data (log defl, log rgdp, short term rate, log credit2gdp, log real property prices, log commodity prices) to use in plutusTvv3...
##          and converts data into an annual time series that begins with year 1
##Input: 
	##countries: vector of country abbreviations
	##data: spreadsheet of current price gdp (<country abbrev>_gdp), gdp deflator (<country abbrev>_defl), private credit (<country abbrev>_credittype),
	##      short term interest rate (<country abbrev>_rate), real residential prices (rp_<country abbrev>)
	##commdata: spreadsheet of commodity index prices (USD) and country exchange rates (USD/LCU)). A World Bank monthly time series. 
		#For J:
		#commdata=read.csv("comm.csv",  na.strings="NA")
		#commdata<-ts(commdata, start=c(1971,1), freq=12). Done before.
	##commIndexWeights =(0,0,0,0) on the indices (iENERGY ,iNONFUEL,iPRECIOUSMET, copper)	
	##freq: frequency of the data (default=4, quarterly)
	##nLags: nLags in VAR 	
	##SecondaryBrk: dates before which sigmas change
	##begin: date all country time series should be set to start. NULL allows time series to start as far back as possible.
	##terminate: date all country time series should be set to end. NULL allows time series to end as late as possible.
	##credittype: BIS credit series ("nbpa", "aha", "ahu", "ana", "apa", "apu", "bpa", "bpu")
	##power of: exponent on credit data. Default is set to one.
##Output:
	##Y: T x nvar dependent variable data matrix. 
	##X: T x nx exogenous variable data matrix. This is the country dummy matrix
	##Tsigbrks: Dates before which sigmas change
	##Tcountrybrk: Dates before which you enter another country's data series. These are "breaks" in the data.
	##countries: List of countries (repeats country name by number of time breaks) (omits countries that didn't have the credit series desired)
	##Start: start dates of country series
	##End: end dates of country series
	##SecondaryBrk: Spits out the Secondary Brk input without the first and last dates
	##lmdseq: Index of the first SecondaryBrk after which the data series begins. For example, if the secondary breaks are 0, 1984, 2006, 2010, a lmdseq=2 means the
	##        data series begins after 0, and before 1984.		
	##lmdendseq: Index of the Secondary Brk immediately before/at which the data series ends	.  With our default settings, a lmdendseq=3 would denote that the data series 
	##           ended on or before 2006.

		
#Prepare the data
###
nCountries0=length(countries)
nCountries=length(countries)

###
#Lists
###
Countrylst=array()
Y<-list()
Xffr<-list()
Startdates=list()
Enddates=list()
dlength0=vector("numeric")
dlength=vector("numeric")
lmdseq=list()
lmdendseq=list()


for (i in 1:nCountries){

		

	iCountry=countries[i]
	a<- paste(iCountry, "gdp", sep="_")
	b<- paste(iCountry, "defl", sep="_")
	
	d<- paste(iCountry, "rate", sep="_")
	f<- paste("rp", iCountry, sep="_")
	
	#exchange rates are USD/LCU
	g<-paste("ex", iCountry, sep="_")

#Credit2GDP data	
	
	credit=NULL
	
		if(!credittype=="nbpa"){

			c<-paste(iCountry, credittype, sep="_")
	
			try(credit <- data[,c], silent=TRUE)
			
			}		
			
		if(credittype=="nbpa"){
			
			g<-paste(iCountry, "apa", sep="_")
			h<-paste(iCountry, "bpa", sep="_")
			
			credit1=NULL
			credit2=NULL
			
			try(credit1<-data[,g], silent=TRUE)
			try(credit2<-data[,h], silent=TRUE)
		
		#knowingfully misnaming credit to credit2gdp
			if (!is.null(credit1) & !is.null(credit2)){		
			credit <- credit1-credit2
			}
			}		
			
#If country has credit series data, concatenate the transformed data	
	if(!is.null(credit)){
		
	credit2gdp=log(credit/data[,a])

#Transform rgdp, strate, defl, rrpi data		
	rgdp <- log(data[,a]/data[,b])
	strate <- 100*data[,d]
	defl <- log(data[,b])
	rrpi<-log(data[,f])
	ffr<- data[, "us_ffr"]
	
	datay <- cbind(defl, rgdp, credit2gdp, rrpi, strate)
	
	
#commodityIndex Data (in USD)
	comm <- log(commdata[,g]*(commIndexWeights[1]*commdata[,"iENERGY"] + commIndexWeights[2]*commdata[,"iNONFUEL"] + commIndexWeights[3]*commdata[,"iPRECIOUSMET"]+commIndexWeights[4]*commdata[,"Copper"]))
	comm <- aggregate(comm, nfreq=4)

#Time set the data					

	datay <- ts(datay, start=c(1955,1), freq=freq)
	
	#Combine with commodities data
	datay <- ts.union(datay, comm)
	
	#Only take observations for which there exists observations for all variables
	datay <-trimts(datay)
	datay <- diff(datay)
	
	varnames=c("Defl", "rGDP", "C2G", "rPP", "IR", "CP")
	dimnames(datay)[[2]]=varnames
	dimnames(datay)[[2]]=varnames
	nVar<-dim(datay)[[2]]
	
	
#Specify start or end date of the data if desired
	if (!is.null(begin)){
	datay = window(datay, start=begin)
						}	
						
	if (!is.null(terminate)){
		datay = window(datay, end=terminate)
						}
						
						
#Generate the dates for country breaks
	dlength0[i]<-dim(datay)[1]	
							
							
#Fill in lmdseq, which gives you the index of the first SecondaryBrk before which the data series begins 
	lmdseq[[i]]=NA
	lmdendseq[[i]]=NA
	nSecondaryBrks=length(SecondaryBrk)
	
	
	for (m in 1:nSecondaryBrks){
		if (is.na(lmdseq[[i]])){
		if (tsp(datay)[1]<SecondaryBrk[m]){lmdseq[[i]]=m-1}
		}
		if (is.na(lmdendseq[[i]])){
		if (tsp(datay)[2]<=SecondaryBrk[m]){lmdendseq[[i]]=m}
		}
		}

	#Generate length of every country-time sub series
	for (j in (lmdseq[[i]]:(lmdendseq[[i]]-1))){
	datay0=NULL
	datay0=window(datay, start=SecondaryBrk[j]+ 1/freq, end=SecondaryBrk[j+1])				
	Y[[(i-1)*(nSecondaryBrks-1)+j]]=datay0
	dlength[[(i-1)*(nSecondaryBrks-1)+j]]=dim(datay0)[1]
	Startdates[[(i-1)*(nSecondaryBrks-1)+j]]=tsp(datay0)[1]
	Enddates[[(i-1)*(nSecondaryBrks-1)+j]]=tsp(datay0)[2]
	Countrylst[[(i-1)*(nSecondaryBrks-1)+j]]=iCountry
	}
	}
	}


###
#Prep Countrylst 
###
	Countrylst=na.omit(Countrylst)

	
##
#Prep dlength0 (used to get Tcountrybrk), and dlength (used to get Tsigbrk)
###

	dlength0=na.omit(dlength0)
	dlength=na.omit(dlength)
	nBrks =length(dlength)

##
#Prep Y Matrix 
###
	
	Y<-do.call("rbind", Y)
	Y=ts(Y, start=1, frequency=1)
	Tlength=dim(Y)[1]


##
#Prep X Matrix--matrix of country dummies
###
	
	X <- matrix(0,Tlength, nCountries0)
	start <- 0
	for (ic in 1:(nCountries)) {
		X[(start + 1:dlength0[ic]), ic] <- 1
		start <- start + dlength0[ic]
		}
	
##
#Prep X Matrix--matrix of country dummies
###	
	#start<-0
	#X2<-matrix(0, Tlength, nCountries0*nVar*(nLags+1))
	
	#for (ic in 1:nCountries0) {
	#		for (iLag in 0:nLags){
	#		X2[start + (nLags+1):dlength0[ic], (ic-1)*nVar*(nLags+1) + nVar*(iLag)+ (1:nVar)] <-Y[start -iLag + (nLags+1):dlength0[ic],]
	#		}
	#	start <- start + dlength0[ic]
	#	}	
	#X<-cbind(X,X2)
		
###
#Prep Tsigbrk
###
	dlength=dlength[1:(nBrks-1)]
	Countrybrk <- cumsum(dlength0)
	Tsigbrk <- cumsum(dlength)
	
	Startdates=unlist(Startdates)
	names(Startdates)=Countrylst
	
	Enddates=unlist(Enddates)
	names(Enddates)=Countrylst
	
	SecondaryBrk=SecondaryBrk[2:(nSecondaryBrks-1)]
	return(list(Y=Y,X=X, Countrybrk = Countrybrk ,Tsigbrk=Tsigbrk, countries=Countrylst, dlength=dlength, dlength0=dlength0, Start=Startdates,End=Enddates, SecondaryBrk=SecondaryBrk, lmdseq=lmdseq, lmdendseq=lmdendseq))
}

################################
#Execute
################################

