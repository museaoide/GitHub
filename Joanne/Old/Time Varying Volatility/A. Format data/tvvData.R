
################################
#Function
################################

tvvData<-function(data, countries, begin=c(1985,4), terminate=c(2008,1), credittype="pc", powerof=1){

	#function: 
		#stacks country data (log defl, log rgdp, short term rate, log credit2gdp, log real, property prices) to use in plutusTvv.., and converts data into an annual time series that begins with year 1
		#ffr is an exogenous variable
	#data: spreadsheet of current price gdp (country abbrev_gdp), gdp deflator (country abbrev_defl), private credit (country abbrev_credittype), rate (country abbrev_rate), real residential prices (rp_country abbrev)
	#countries: vector of country abbreviation.
	#credittype: BIS credit series ("nbpa", "aha", "ahu", "ana", "apa", "apu", "bpa", "bpu")
	#power of: change if you want to estimate a nonlinear model for credit2gdp. Default is 1
	#output: list
		#Y: T x nvar dependent variable data matrix. 
		#X: T x nx exogenous variable data matrix. This is the country dummy matrix
		#Tsigbrks: "structural breaks", which for us is the last observation date for a country
		#countries: countries (omits countries that didn't have the credit series desired)
		#dlength: length of series by country
		#Start: start date of a country's series
		#End: enddate of a country's series

##
#Prepare the data
##
nCountries=length(countries)
Countrylst=array()
Y<-list()
Xffr<-list()
Startdates=list()
Enddates=list()
dlength=vector("numeric")


for (i in 1:nCountries){

	iCountry=countries[i]
	a<- paste(iCountry, "gdp", sep="_")
	b<- paste(iCountry, "defl", sep="_")
	
	d<- paste(iCountry, "rate", sep="_")
	f<- paste("rp", iCountry, sep="_")
	
	
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
			
		
	
	rgdp <- log(data[,a]/data[,b]/100)
	strate <- 100*data[,d]
	defl <- log(data[,b])
	rrpi<-log(data[,f])
	ffr<- (data[, "us_ffr"])
	
	if(!is.null(credit2gdp)){
		
	credit2gdp=log(credit2gdp/data[,a])
	datay <- cbind(defl, rgdp, credit2gdp, rrpi, strate, ffr)
	datay <- ts(datay, start=c(1955,1), freq=4)
	datay <-trimts(datay)
			
		#seasonal=decompose(datay[, "credit2gdp"])
		#deseasonal= seasonal$x - seasonal$seasonal
	
		#datay[,"credit2gdp"]= log(deseasonal/datay[,"rgdp"])


	if (!is.null(begin) && !is.null(terminate)){
	datay = window(datay, start=begin, end=terminate)
	}
	
	start=tsp(datay)[1]
	end=tsp(datay)[2]
	datay <- ts(datay, start= 1, freq=1)
	datax <- datay[,6]
	datay <- datay[, 1:5]
	
	Y[[i]]<- datay
	Xffr[[i]]<-datax
	Startdates[[iCountry]]= start
	Enddates[[iCountry]]=end
	
	dlength[i]<-dim(datay)[1]
	Countrylst[i]=iCountry
	}
	
	}

##
#Prep Countrylst 
###
Countrylst=na.omit(Countrylst)

	
##
#Prep Dlength 
###
dlength=na.omit(dlength)
nCountries=length(dlength)

##
#Prep Y Matrix
###

Xffr=ts(unlist(Xffr), start=1, frequency=1)

Y<-do.call("rbind", Y)
Y=ts(Y, start=1, frequency=1)
Tlength=dim(Y)[1]

##
#Prep X Matrix
###
	
 X <- matrix(0,Tlength, nCountries)
 start <- 0
for (ic in 1:nCountries) {
	X[(start + 1:dlength[ic]), ic] <- 1
	start <- start + dlength[ic]
	}
	
 X=cbind(X, Xffr)

##
#Prep Tsigbrk
###
dlength=dlength[1:(nCountries-1)]
Tsigbrk <- cumsum(dlength)

return(list(Y=Y,X=X, Tsigbrk=Tsigbrk, countries=Countrylst, dlength=dlength, Start=unlist(Startdates), End=unlist(Enddates)))
}

################################
#Execute
################################
