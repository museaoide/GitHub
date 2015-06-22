summApos<-function(countries, resultslst){
nCountries=length(countries)
nVar=6
summat=matrix(0, nVar, nVar)

for (i in 1:nCountries){
	iCountry=countries[i]
	countrymatA=resultslst[[i]]$A
		for (j in 1:nVar){
		for (k in 1:nVar){
			if(countrymatA[j,k]>0){
			summat[j,k]=summat[j,k]+1	
			}
		}
	} 
}

return(summat)
}

summAneg<-function(countries, resultslst){
nCountries=length(countries)
nVar=6
summat=matrix(0, nVar, nVar)

for (i in 1:nCountries){
	iCountry=countries[i]
	countrymatA=resultslst[[i]]$A
		for (j in 1:nVar){
		for (k in 1:nVar){
			if(countrymatA[j,k]<0){
			summat[j,k]=summat[j,k]+1	
			}
		}
	} 
}

return(summat)
}

countries=c("aus", "be", "can", "che", "dan", "fr", "it", "jpn", "nl", "no", "pt", "se", "uk", "us", "za")

resultslst=list()
resultslst[[1]]<-results6e2_can
resultslst[[2]]<-results6e2_aus
resultslst[[3]]<-results6e2_be
resultslst[[4]]<-results6e2_che
resultslst[[5]]<-results6e2_dan
resultslst[[6]]<-results6e2_fr
resultslst[[7]]<-results6e2_it
resultslst[[8]]<-results6e2_jpn
resultslst[[9]]<-results6e2_nl
resultslst[[10]]<-results6e2_no
resultslst[[11]]<-results6e2_pt
resultslst[[12]]<-results6e2_se
resultslst[[13]]<-results6e2_uk
resultslst[[14]]<-results6e2_us
resultslst[[15]]<-results6e2_za
resultslst[[16]]<- results6e2_be
resultslst[[17]]<-results6e2_fr
resultslst[[18]]<-results6e2_it
resultslst[[19]]<-results6e2_nl
resultslst[[20]]<-results6e2_pt


lcA0=matrix(TRUE, 6,6)
diag(lcA0)=FALSE
seedA=NULL
nCountries=length(unique(listData6e2_14$countries))
for (i in 1:nCountries){
seedA=c(seedx, resultslst[[i]]$A[lcA0])
}


summApos(countries, resultslst)
summAneg(countries, resultslst)
