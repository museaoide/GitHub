


##################################FUNCTION: RUN VARS, IRs
#12/1/2004

rfvar7 <- function(results, nSteps=60){
	#function: runs reduced form VARs on residuals from the international SVAR model
	#results: results from internatl SVAR model
	
	u=results$vout$var$u
	countries=results$opt$listData$countries
	nCountries=length(countries)
	regimes=results$regimes	
	voutlst<-list()
	
for (i in 1:(nCountries-1)){
	iCountry=countries[i]
	iRegime= which(iCountry==countries, arr.ind=TRUE)
	iRange=(regimes[iRegime]+1):(regimes[iRegime+1])
	iU=u[iRange,]
	dimnames(iU)=dimnames(results1$opt$listData$Y)
	
	#reduced form VAR		
		vout <-rfvar3(ydata=iU, lags=5, const=TRUE)
		voutlst[[iCountry]]<-vout	
}


	i=nCountries
	iCountry=countries[i]
	iRegime= which(iCountry==countries, arr.ind=TRUE)
	iU=u[(regimes[iRegime]+1): dim(iU)[1],]
	voutlst[[iCountry]]<-vout	


	return(list(voutlst=voutlst, countries=countries))
	}



tvvData_noffr2 <- function(rfvar7results){
	#function: runs the international model on the reduced form VARs

countries=rfvar7results$countries	
nCountries=length(countries)
dlength<-list()
Y<-list()

for (i in 1:nCountries){
	iCountry=countries[i]
	datay= rfvar7results$voutlst[[iCountry]]$u
	Y[[iCountry]]=datay
	dlength[i]=dim(datay)[1]
	}
	
	Y <- do.call("rbind", Y)
	Y=ts(Y, start=1, frequency=1)
	Tlength=dim(Y)[1]
 	dlength=unlist(dlength)
 	
	X <- matrix(0,Tlength, nCountries)
		start <- 0
		for (ic in 1:nCountries) {
			X[(start + 1:dlength[ic]), ic] <- 1
			start <- start + dlength[ic]
			}
	
	
	dlength=dlength[1:(nCountries-1)]
	Tsigbrk=cumsum(dlength)
	
	return(list(Y=Y, X=X, Tsigbrk=Tsigbrk, countries=countries))
	}

#Jointoptimization

#We need to think of a criteria for the "termination" of the search
#You can do it with a joint liklihood, or you could see if the coefficients are near zero.

#run reduced form	
rvar7results_2=rfvar7(results1c)
results1c_2=plutusTvv(tvvData_noffr2(rvar7results_2), H1=10^-6, H2=10^-8)
save(results1c_2, file="results1c.RData")





