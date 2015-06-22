
listData2 <- function(rfvar7results){
	#function: runs the international model on the reduced form VARs

countries=rfvar7results$countries	
nCountries=length(countries)
Y<-list()

for (iCountry in countries){
	datay= rfvar7results$voutlst[[iCountry]]$u
	
	Y[[i]]=datay
	dlength[i]=dim(datay)
	}
	
	Y <- do.call("rbind", Y)
	Y=ts(Y, start=1, frequency=1)
	Tlength=dim(Y)[1]
 
	X <- matrix(0,Tlength, nCountries)
		start <- 0
		for (ic in 1:(nCountries) {
			X[(start + 1:dlength[ic]), ic] <- 1
			start <- start + dlength[ic]
			}
	
	dlength=dlength[1:(nCountries-1)]
	Tsigbrk =cumsum(dlength)
	
	return(list(Y=Y, X=X, Tsigbrk=Tsigbrk, countries=countries)

