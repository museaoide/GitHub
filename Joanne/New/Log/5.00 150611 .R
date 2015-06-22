
*Medium size set (Joanne)
listData_15610a=tvvData_noffr4(countries=c("can", "che", "uk") ,data,commdata, SecondaryBrk=c(0, 2007.75, 2015))
results_15610a=plutusTvv5(listData_15610a,nit=700,xmnprior=list(tight=5,decay=.5))

esresults_15610a= plutusError(endct=50)


*Larger size set. Sims.
##Tight=10
listData_15610c=tvvData_noffr4(countries=c("can", "che", "uk", "jpn", "aus", "za") ,data,commdata, SecondaryBrk=c(0, 2015))
results_15610c=plutusTvv5(listData_15610c,nit=700,xmnprior=list(tight=10,decay=.5), seedx=x1)
##Tight=20
listData_15610d=tvvData_noffr4(countries=c("can", "che", "uk", "jpn", "aus", "za") ,data,commdata, SecondaryBrk=c(0, 2015))
results_15610d=plutusTvv5(listData_15610c,nit=700,xmnprior=list(tight=20,decay=.5), seedx=x1)
##Tight=2
listData_15610e=tvvData_noffr4(countries=c("can", "che", "uk", "jpn", "aus", "za") ,data,commdata, SecondaryBrk=c(0, 2015))
results_15610e=plutusTvv5(listData_15610c,nit=700,xmnprior=list(tight=2,decay=.5), seedx=x1)

*Decreased tightness of the mnprior, which was set at tight=3.
listData_15610f=tvvData_noffr4(countries=c("can", "che", "uk", "jpn", "aus", "za") ,data,commdata, SecondaryBrk=c(0, 2015))
results_15610f=plutusTvv5(listData_15610c,nit=700,mnprior=list(tight=1, decay=.5), xmnprior=list(tight=2,decay=.5), seedx=x1)

*Decreased tightness of mn prior and increased tightness of xmnprior for comparison with *f (1:35)
results_15610g=plutusTvv5(listData_15610c,nit=700,mnprior=list(tight=1, decay=.5), xmnprior=list(tight=10,decay=.5), seedx=x1)

results_15610c= plutusError(endct=50)

*
results_15610h=plutusTvv5(listData_15610c,nit=700,xmnprior=list(tight=5,decay=.5), seedx=x1)


--
##Allow for country specific impulse response

##The first dimension is the equation number
##By = nvar x nvar x lags

##Create By2 matrix
Bx = nvar x nx
By2= Bx[nvar, 1:(nc*nv*lags)]
dim(By2)= c(nv, nC , vars, lags)

for (ic in 1: nc){
	iBy= iBBy + By2[ ,ic, , ]
}


--

Things to write up: 
1. Explain change in units (make sure prior variances are of a correct scale, thso that implementation of our beliefs are correct i magnitude)
Things to do
1.) Iterate over a number of lambda candidates
2.) 