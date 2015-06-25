
##
#Completed
##

*NO MP SHOCK. BAD.
listData_15610c=tvvData_noffr4(countries=c("can", "che", "uk", "jpn", "aus", "za") ,data,commdata, SecondaryBrk=c(0, 2015))
results_15610d=plutusTvv5(listData_15610c,nit=700,xmnprior=list(tight=20,decay=.5), seedx=x1)


*NO MP SHOCK. Looking OK
results_15610e=plutusTvv5(listData_15610c,nit=700,xmnprior=list(tight=2,decay=.5), seedx=x1)

*NO MP SHOCK. Medium Bad.
results_15610h=plutusTvv5(listData_15610c,nit=700,xmnprior=list(tight=5,decay=.5), seedx=x1)

*NO MP Shock. Medium Bad.
results_15610i=plutusTvv5(listData_15610c,nit=700,xmnprior=list(tight=7,decay=dx=x1)

*NO MP SHOCK. BAD.
results_15610f=plutusTvv5(listData_15610c,nit=700,mnprior=list(tight=1, decay=.5), xmnprior=list(tight=2,decay=.5), seedx=x1)

##
#RUNNING at 1:37 A.M.
##

*NEXT. Try tigher Minnesota Prior, and looser prior on the prior over the country ys.
results_15610j=plutusTvv5(listData_15610c,nit=700,mnprior=list(tight=4, decay=.5), xmnprior=list(tight=1,decay=.5), seedx=x1)


*NEXT. Try tigher Minnesota Prior, and tighter prior on the prior over the country ys.
results_15610k=plutusTvv5(listData_15610c,nit=700,mnprior=list(tight=4, decay=.5), xmnprior=list(tight=2,decay=.5), seedx=x1)


listData_15610c=tvvData_noffr4(countries=c("can", "che", "uk", "jpn", "aus", "za") ,data,commdata, SecondaryBrk=c(0, 2015))


results_15610l=plutusTvv5(listData_15610c,nit=700,mnprior=list(tight=5, decay=.5), xmnprior=list(tight=2,decay=.5), seedx=x1)
results_15610m=plutusTvv5(listData_15610c,nit=700,mnprior=list(tight=3, decay=.5), xmnprior=list(tight=1,decay=.5), seedx=x1)
results_15610n=plutusTvv5(listData_15610c,nit=700,mnprior=list(tight=5, decay=.5), xmnprior=list(tight=5,decay=.5), seedx=x1)


##
#TO DO
##
Re-run, and do mcmc draws to have results to compare to. Something is wrong in the old code, and I have to figure out what.


##
#Write up
##

Things to write up: 
1. Explain change in units (make sure prior variances are of a correct scale, thso that implementation of our beliefs are correct i magnitude)
Things to do
1.) Iterate over a number of lambda candidates
2.) 
