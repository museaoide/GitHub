results
##################################Function


IPunitrootdiscard<- function(results, countries, credittype="pc", powerof=1, nSteps=40, nDraws=1000, respvar="rgdp", shockvar="credit2gdp", shocksize=1, filename="impulse_plot", errorbands=TRUE, cinterval=list( c(.16, .84), c(.025, .975)), maxdrawpct=.5, unitroot=1.0025, scale=TRUE){

#function: modification of IPunitrootdiscard that discards posterior draws with unit roots
#data: output from plutus4vars/plutus5vars
#countries: vector of country abbreviations
#credittype: bis credit series
#powerof: power on credit-to-gdp 
#nSteps: number of steps for IR response
#nDraws: number of trials for posterior draw
#resp var: response variable for IR response
#shock var: shock variable for IR response
#shock size: size of shock
#cintervals: intervals for the error bands. Should be ordered from smallest to largest. 
#maxdrawpct, unitroot are inputs for postdrawDiscard
#scale=TRUE if y-axis scaled to 1 std dev shock

#define plot area
nCountries=length(countries)
arr = c(ceiling(nCountries/6), 6)
par(mfrow=arr, col.lab="black", col.main="black", oma=c(1,3,3,2), mar=c(1, 2, 1, 2), tcl = -.1, mgp=c(0,0,0))
axisline=rep(0, (nSteps + 2))

#pull data
for (iCountry in countries){
resp <- results[[2]][[iCountry]]	
vout<- results[[1]][[iCountry]]
if (errorbands==TRUE){	

###--Different from impulseplot2 here	
	postdraws<-postdrawDiscard(vout, nDraws=nDraws, returnfrac=FALSE,maxdrawpct=maxdrawpct, unitroot=unitroot)
###-----
	
#Posterior simulation IR
	nDraws=dim(postdraws$smat)[3]
	irdraws=array(NA, c(dim(resp), nDraws))
	dimnames(irdraws) <- list(colnames(data[[2]][[iCountry]]), colnames(data[[2]][[iCountry]]), NULL, NULL)
	for (iDraw in (1:nDraws)){
		drawReg=list(By=postdraws$By[,,,iDraw], Bx=postdraws$Bx[,iDraw])
		irdraws[,,,iDraw]=impulsdtrf(drawReg, smat=t(postdraws$smat[,,iDraw]), nstep=dim(irdraws)[3])
	}
#Draw error bands	
	Resptrials=shocksize*irdraws[respvar,,,]
	Shocktrials=Resptrials[shockvar, , ]
	
	Respbounds<-list()	
	nCintervals=length(cinterval)
	gRange = 1:(nSteps)	
for (iCinterval in 1: nCintervals){	
	iShockbounds<-t(apply(Shocktrials, 1, quantile, probs=c(cinterval[[iCinterval]][1], cinterval[[iCinterval]][2])))
	iRespbounds<-array(NA, c(nSteps, 2, 1))
	iRespbounds[,,1] = iShockbounds[1:(nSteps),]
	Respbounds[[iCinterval]]<-iRespbounds
	}
	}
else {
	Respbounds=NA
	}			
#Draw impulse response
Respseries=shocksize*resp[respvar, shockvar,]	
#Set up to find graph range	
yMax = max(Respseries, unlist(Respbounds), 0, na.rm=TRUE) + .0001 * abs(max(Respseries, unlist(Respbounds),na.rm=TRUE))
yMin =  min(Respseries, unlist(Respbounds), 0,na.rm=TRUE) - .0001 * abs(min(Respseries, unlist(Respbounds),na.rm=TRUE))	
yMax=signif(yMax, 3)
yMin=signif(yMin,3)
if (yMin<0 & yMax>0){
	if (abs(yMin)>abs(yMax)){
		yMax=abs(yMin)}
	else{yMin=-yMax}}
	
if (scale==TRUE){
iUstd=vout$u
Omega = crossprod(iUstd) / dim(iUstd)[1]
sdev = sqrt(Omega[shockvar, shockvar])	
		
yMin=yMin/sdev
yMax=yMax/sdev
}
	


#plot series	
title=iCountry
plot(1:nSteps, Respseries[1:nSteps], ylim=c(yMin, yMax), type="l", main=title, lwd=1, xlab='', ylab='', yaxt='n', xaxt='n')

if (errorbands==TRUE){
upper<-list()
lower <- list()		
for (iCinterval in 1: nCintervals)	{
iupper<-Respbounds[[iCinterval]][,2,1]
upper[[iCinterval]]<-iupper
ilower<-Respbounds[[iCinterval]][,1,1]	
lower[[iCinterval]]<- ilower
}
for (iCinterval in 1: nCintervals)	{	
if (!is.null(upper[[iCinterval]])) {
			polygon(c(gRange, rev(gRange)), c(lower[[iCinterval]][gRange], rev(upper[[iCinterval]][gRange])),  
			col=rgb(0, 0, 1/iCinterval,0.3))}}}

			
lines(axisline, type='l', lwd=1)
axis(side=2, at=c(yMin, 0, yMax), labels=TRUE, las=2,cex.axis=.7, tck=.03)
	}

bigtitle = paste('Impulse response of', respvar, 'to', shockvar, 'shock over', as.character(nSteps), 'quarters', sep = ' ')

title(bigtitle, outer = TRUE, cex = 1.2)

dev.copy2pdf(file = filename)}

##################################EXECUTE	


IPunitrootdiscard(plutus5vars(data, countriesnur1.25), countriesnur1.25, respvar="rgdp", shockvar="credit2gdp", errorbands=TRUE, cinterval=list( c(.16, .84), c(.025, .975)), filename="IR_rgdp_c2gshock_68_95", unitroot=1.25)

IPunitrootdiscard(plutus5vars(data, countriesnur1.025), countriesnur1.025, respvar="rgdp", shockvar="credit2gdp", errorbands=TRUE, cinterval=list( c(.16, .84), c(.025, .975)), filename="IR_rgdp_c2gshock_68_95", unitroot=1.025)

IPunitrootdiscard(plutus5vars(data, countriesnur1.01), countriesnur1.01, respvar="rgdp", shockvar="credit2gdp", errorbands=TRUE, cinterval=list( c(.16, .84), c(.025, .975)), filename="IR_rgdp_c2gshock_68_95", unitroot=1.01, scale=FALSE)


#######
With Error bands
######
###
Resp of rgdp
###

IPunitrootdiscard(results=red, countries=countriesffr, respvar="rgdp", shockvar="credit2gdp", errorbands=TRUE, cinterval=list( c(.16, .84), c(.025, .975)), filename="IR_rgdp_c2gshock_68_95")
IPunitrootdiscard(results=red, countries=countriesffr, respvar="rgdp", shockvar="defl", errorbands=TRUE, cinterval=list( c(.16, .84), c(.025, .975)), filename="IR_rgdp_defl_68_95")
IPunitrootdiscard(results=red, countries=countriesffr, respvar="rgdp", shockvar="ffr", errorbands=TRUE, cinterval=list( c(.16, .84), c(.025, .975)), filename="IR_rgdp_ffr_68_95")


###
Resp of credit2gdp
###


IPunitrootdiscard(plutus5vars(data, countries), countries, respvar="credit2gdp", shockvar="rgdp", errorbands=TRUE, cinterval=list( c(.16, .84), c(.025, .975)), filename="IR_c2g_rgdp_68_95")
IPunitrootdiscard(plutus5vars(data, countries), countries, respvar="credit2gdp", shockvar="defl", errorbands=TRUE, cinterval=list( c(.16, .84), c(.025, .975)), filename="IR_c2g_defl_68_95")
IPunitrootdiscard(plutus5vars(data, countries), countries, respvar="credit2gdp", shockvar="strate", errorbands=TRUE, cinterval=list( c(.16, .84), c(.025, .975)), filename="IR_c2g_strate_68_95")
IPunitrootdiscard(plutus5vars(data, countries), countries, respvar="credit2gdp", shockvar="ffr", errorbands=TRUE, cinterval=list( c(.16, .84), c(.025, .975)), filename="IR_c2g_ffr_68_95")

###
Resp of defl
###


IPunitrootdiscard(plutus5vars(data, countries), countries, respvar="defl", shockvar="rgdp", errorbands=TRUE, cinterval=list( c(.16, .84), c(.025, .975)), filename="IR_defl_rgdp_68_95")
IPunitrootdiscard(plutus5vars(data, countries), countries, respvar="defl", shockvar="credit2gdp", errorbands=TRUE, cinterval=list( c(.16, .84), c(.025, .975)), filename="IR_defl_credit2gdpe_68_95")
IPunitrootdiscard(plutus5vars(data, countries), countries, respvar="defl", shockvar="strate", errorbands=TRUE, cinterval=list( c(.16, .84), c(.025, .975)), filename="IR_defl_strate_68_95")
IPunitrootdiscard(plutus5vars(data, countries), countries, respvar="defl", shockvar="ffr", errorbands=TRUE, cinterval=list( c(.16, .84), c(.025, .975)), filename="IR_defl_ffr_68_95")

###
Resp of strate
###


IPunitrootdiscard(plutus5vars(data, countries), countries, respvar="strate", shockvar="rgdp", errorbands=TRUE, cinterval=list( c(.16, .84), c(.025, .975)), filename="IR_strate_rgdp_68_95")
IPunitrootdiscard(plutus5vars(data, countries), countries, respvar="strate", shockvar="credit2gdp", errorbands=TRUE, cinterval=list( c(.16, .84), c(.025, .975)), filename="IR_strate_credit2gdp_68_95")
IPunitrootdiscard(plutus5vars(data, countries), countries, respvar="strate", shockvar="defl", errorbands=TRUE, cinterval=list( c(.16, .84), c(.025, .975)), filename="IR_strate_defl_68_95")
IPunitrootdiscard(plutus5vars(data, countries), countries, respvar="strate", shockvar="ffr", errorbands=TRUE, cinterval=list( c(.16, .84), c(.025, .975)), filename="IR_strate_ffr_68_95")

IPunitrootdiscard(plutus5vars(data, countries), countries, respvar="ffr", shockvar="rgdp", errorbands=TRUE, cinterval=list( c(.16, .84), c(.025, .975)), filename="IR_ffr_rgdp_68_95")
IPunitrootdiscard(plutus5vars(data, countries), countries, respvar="ffr", shockvar="credit2gdp", errorbands= TRUE, cinterval=list( c(.16, .84), c(.025, .975)), filename="IR_ffr_credit2gdp_68_95")
IPunitrootdiscard(plutus5vars(data, countries), countries, respvar="ffr", shockvar="defl", errorbands= TRUE, cinterval=list( c(.16, .84), c(.025, .975)), filename="IR_ffr_defl_68_95")
IPunitrootdiscard(plutus5vars(data, countries), countries, respvar="ffr", shockvar="strate", errorbands= TRUE, cinterval=list( c(.16, .84), c(.025, .975)), filename="IR_ffr_strate_68_95")


######
Without Error Bands
######

###
Resp of rgdp
###

IPunitrootdiscard(plutus5vars(data, countries), countries, respvar="rgdp", shockvar="credit2gdp", errorbands=FALSE, cinterval=list( c(.16, .84), c(.025, .975)), filename="IR_rgdp_c2gshock_68_95_nobands")
IPunitrootdiscard(plutus5vars(data, countries), countries, respvar="rgdp", shockvar="defl", errorbands= FALSE, cinterval=list( c(.16, .84), c(.025, .975)), filename="IR_rgdp_defl_68_95_nobands")
IPunitrootdiscard(plutus5vars(data, countries), countries, respvar="rgdp", shockvar="strate", errorbands= FALSE, cinterval=list( c(.16, .84), c(.025, .975)), filename="IR_rgdp_strate_68_95_nobands")
IPunitrootdiscard(plutus5vars(data, countries), countries, respvar="rgdp", shockvar="ffr", errorbands= FALSE, cinterval=list( c(.16, .84), c(.025, .975)), filename="IR_rgdp_ffr_68_95_nobands")


###
Resp of credit2gdp
###


IPunitrootdiscard(plutus5vars(data, countries), countries, respvar="credit2gdp", shockvar="rgdp", errorbands= FALSE, cinterval=list( c(.16, .84), c(.025, .975)), filename="IR_c2g_rgdp_68_95_nobands")
IPunitrootdiscard(plutus5vars(data, countries), countries, respvar="credit2gdp", shockvar="defl", errorbands= FALSE, cinterval=list( c(.16, .84), c(.025, .975)), filename="IR_c2g_defl_68_95_nobands")
IPunitrootdiscard(plutus5vars(data, countries), countries, respvar="credit2gdp", shockvar="ffr", errorbands= FALSE, cinterval=list( c(.16, .84), c(.025, .975)), filename="IR_c2g_ffr_68_95_nobands")

###
Resp of defl
###


IPunitrootdiscard(plutus5vars(data, countries), countries, respvar="defl", shockvar="rgdp", errorbands= FALSE, cinterval=list( c(.16, .84), c(.025, .975)), filename="IR_defl_rgdp_68_95_nobands")
IPunitrootdiscard(plutus5vars(data, countries), countries, respvar="defl", shockvar="credit2gdp", errorbands= FALSE, cinterval=list( c(.16, .84), c(.025, .975)), filename="IR_defl_credit2gdpe_68_95_nobands")
IPunitrootdiscard(plutus5vars(data, countries), countries, respvar="defl", shockvar="strate", errorbands= FALSE, cinterval=list( c(.16, .84), c(.025, .975)), filename="IR_defl_strate_68_95_nobands")
IPunitrootdiscard(plutus5vars(data, countries), countries, respvar="defl", shockvar="ffr", errorbands= FALSE, cinterval=list( c(.16, .84), c(.025, .975)), filename="IR_defl_ffr_68_95_nobands")

###
Resp of strate
###


IPunitrootdiscard(plutus5vars(data, countries), countries, respvar="strate", shockvar="rgdp", errorbands=FALSE, cinterval=list( c(.16, .84), c(.025, .975)), filename="IR_strate_rgdp_68_95_nobands")
IPunitrootdiscard(plutus5vars(data, countries), countries, respvar="strate", shockvar="credit2gdp", errorbands= FALSE, cinterval=list( c(.16, .84), c(.025, .975)), filename="IR_strate_credit2gdp_68_95_nobands")
IPunitrootdiscard(plutus5vars(data, countries), countries, respvar="strate", shockvar="defl", errorbands=FALSE, cinterval=list( c(.16, .84), c(.025, .975)), filename="IR_strate_defl_68_95_nobands")
IPunitrootdiscard(plutus5vars(data, countries), countries, respvar="strate", shockvar="ffr", errorbands=FALSE, cinterval=list( c(.16, .84), c(.025, .975)), filename="IR_strate_ffr_68_95_nobands")


###
Resp of ffr
###


IPunitrootdiscard(plutus5vars(data, countries), countries, respvar="ffr", shockvar="rgdp", errorbands=FALSE, cinterval=list( c(.16, .84), c(.025, .975)), filename="IR_ffr_rgdp_68_95_nobands")
IPunitrootdiscard(plutus5vars(data, countries), countries, respvar="ffr", shockvar="credit2gdp", errorbands= FALSE, cinterval=list( c(.16, .84), c(.025, .975)), filename="IR_ffr_credit2gdp_68_95_nobands")
IPunitrootdiscard(plutus5vars(data, countries), countries, respvar="ffr", shockvar="defl", errorbands=TRUE, cinterval=list( c(.16, .84), c(.025, .975)), filename="IR_ffr_defl_68_95_nobands")
IPunitrootdiscard(plutus5vars(data, countries), countries, respvar="ffr", shockvar="strate", errorbands=TRUE, cinterval=list( c(.16, .84), c(.025, .975)), filename="IR_ffr_strate_68_95_nobands")

