################################
#Function
################################

tvvData_noffr2<-function(data, countries, SecondaryBrk=c(0, 1984, 2006, 2010, 2014), begin=NULL, terminate=NULL, credittype="apa", powerof=1){
	##Function: stacks country data (log defl, log rgdp, short term rate, log credit2gdp, log real property prices) to use in plutusTvv3..,
	##and converts data into an annual time series that begins with year 1
		##Data: spreadsheet of current price gdp (country abbrev_gdp), gdp deflator (country abbrev_defl), private credit (country abbrev_credittype),
		##rate (country abbrev_rate), real residential prices (rp_country abbrev)
		##Countries: vector of country abbreviation.
		##Credittype: BIS credit series ("nbpa", "aha", "ahu", "ana", "apa", "apu", "bpa", "bpu")
		##Power of: change if you want to estimate a nonlinear model for credit2gdp. Default is set to one.
	##Output: list
		##Y: T x nvar dependent variable data matrix. 
		##X: T x nx exogenous variable data matrix. This is the country dummy matrix
		##Tsigbrks: Dates on which sigmas changes
		##Tcountrybrk: Dates before which you enter another country's data series. These are "breaks" in the data.
		##countries: List of countries (omits countries that didn't have the credit series desired)
		##Start: start date of a country's series
		##End: enddate of a country's series
		##SecondaryBrk: Spits out the Secondary Brk input
		##lmdseq: list of
		
#Prepare the data
###
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


for (i in 1:nCountries){

	iCountry=countries[i]
	a<- paste(iCountry, "gdp", sep="_")
	b<- paste(iCountry, "defl", sep="_")
	
	d<- paste(iCountry, "rate", sep="_")
	f<- paste("rp", iCountry, sep="_")

#Credit2GDP data	
	
	credit2gdp=NULL
	
		if(!credittype=="nbpa"){

			c<-paste(iCountry, credittype, sep="_")
	
			try(credit2gdp <- data[,c], silent=TRUE)
			
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
			credit2gdp<- credit1-credit2
			}
			}			
			
#Transform rgdp, strate, defl, rrpi data		
	
	rgdp <- log(data[,a]/data[,b]/100)
	strate <- 100*data[,d]
	defl <- log(data[,b])
	rrpi<-log(data[,f])
	ffr<- data[, "us_ffr"]

#If country has credit series data, concatenate the transformed data	
	if(!is.null(credit2gdp)){
		
	credit2gdp=log(credit2gdp/data[,a])
	datay <- cbind(defl, rgdp, credit2gdp, rrpi, strate)
	
	#Time set the data
	datay <- ts(datay, start=c(1955,1), freq=4)
	#Only take obs for which there exists obs for all variables
	datay <-trimts(datay)
		
#If desired, truncate the data 
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
nSecondaryBrks=length(SecondaryBrk)

for (m in 1:nSecondaryBrks){
	if (is.na(lmdseq[[i]])){
	if (tsp(datay)[1]<SecondaryBrk[m]){lmdseq[[i]]=m-1}
	}}
	
#Generate length of every country-time sub series
for (j in lmdseq[[i]]:(nSecondaryBrks-1)){
datay0=NULL
datay0=window(datay, start=SecondaryBrk[j]+.25, end=SecondaryBrk[j+1])				
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
nCountries=length(dlength)

##
#Prep Y Matrix
###

Y<-do.call("rbind", Y)
Y=ts(Y, start=1, frequency=1)
Tlength=dim(Y)[1]

##
#Prep X Matrix--matrix of country dummies
###
	
 X <- matrix(0,Tlength, nCountries)
 start <- 0
for (ic in 1:(nCountries)) {
	X[(start + 1:dlength[ic]), ic] <- 1
	start <- start + dlength[ic]
	}
	
###
#Prep Tsigbrk
###
dlength=dlength[1:(nCountries-1)]
Countrybrk <- cumsum(dlength0)
Tsigbrk <- cumsum(dlength)

Startdates=unlist(Startdates)
names(Startdates)=Countrylst

Enddates=unlist(Enddates)
names(Enddates)=Countrylst

SecondaryBrk=SecondaryBrk[-nSecondaryBrks]
SecondaryBrk=SecondaryBrk[-1]

return(list(Y=Y,X=X, Countrybrk = Countrybrk ,Tsigbrk=Tsigbrk, countries=Countrylst, dlength=dlength, dlength0=dlength0, Start=Startdates,End=Enddates, SecondaryBrk=SecondaryBrk, lmdseq=lmdseq))
}

################################
#Execute
################################

