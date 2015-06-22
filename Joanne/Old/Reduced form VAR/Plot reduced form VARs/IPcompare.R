
##################################FUNCTION: PLOT IRS_pc

IPcompare<- function(data1, data2, data2col="red", countries, credittype="pc", powerof=1, nSteps=40, nTrials=1000, respvar="rgdp", shockvar="credit2gdp", shocksize=1, filename="impulse_plot", cinterval=list( c(.16, .84), c(.025, .975))){

#function: modification of impulsePlot that plots *two* impulse responses based on ouput from plutus4vars and or plutus5vars for each country together on the same graph.
#data1: output from plutus4vars/plutus5vars
#data2: output from plutus4vars/plutus5vars. Output structure must match data1.
#countries: vector of country abbreviations
#credittype: bis credit series
#powerof: power on credit-to-gdp 
#nSteps: number of steps for IR response
#nTrials: number of trials for posterior draw
#resp var: response variable for IR response (must be common to both data1 and data1)
#shock var: shock variable for IR response (must be common to both data1 and data1)
#shock size: size of shock
#cintervals: intervals for the error bands. Should be ordered from smallest to largest. 


#define plot area
nCountries=length(countries)
arr = c(ceiling(nCountries/6), 6)
par (mfrow=arr, col.lab="black", col.main="black", oma=c(1,3,3,2), mar=c(1, 2, 1, 2), tcl = -.1, mgp=c(0,0,0))
axisline=rep(0, (nSteps + 2))

#pull data
for (iCountry in countries){
resp1 <- data1[[3]][[iCountry]]	
resp2 <- data2[[3]][[iCountry]]	
#Draw impulse response
Respseries1=shocksize*resp1[respvar, shockvar,]	
Respseries2=shocksize*resp2[respvar, shockvar,]	
#Set up to find graph range	
yMax = max(Respseries1, Respseries2, 0, na.rm=TRUE) + .0001 * abs(max(Respseries1, Respseries2,0, na.rm=TRUE))
yMin =  min(Respseries1, Respseries2, 0, na.rm=TRUE) - .0001 * abs(min(Respseries1, Respseries2,0, na.rm=TRUE))	
yMax=signif(yMax, 3)
yMin=signif(yMin,3)
if (yMin<0 & yMax>0){
	if (abs(yMin)>abs(yMax)){
		yMax=abs(yMin)}
	else{yMin=-yMax}}

#plot series	
title=iCountry
plot(1:nSteps, Respseries1[1:nSteps], ylim=c(yMin, yMax), type="l", main=title, lwd=1, xlab='', ylab='', yaxt='n', xaxt='n')
lines(1:nSteps, Respseries2[1:nSteps], type="l", col= data2col)			
lines(axisline, type='l', lwd=1)
axis(side=2, at=c(yMin, 0, yMax), labels=TRUE, las=2,cex.axis=.7, tck=.03)
	}

bigtitle = paste('Impulse response of', respvar, 'to', shockvar, 'shock over', as.character(nSteps), 'quarters', sep = ' ')

title(bigtitle, outer = TRUE, cex = 1.2)

dev.copy2pdf(file = filename)}

##################################EXECUTE	

###
Resp of rgdp
###
IPcompare(plutus4vars(data, countries), plutus5vars(data, countries), countries=countries, respvar="rgdp", shockvar="credit2gdp", filename="IR_rgdp_credit2gdp_compare")
IPcompare(plutus4vars(data, countries), plutus5vars(data, countries), countries=countries, respvar="rgdp", shockvar="strate", filename="IR_rgdp_strate_compare")
IPcompare(plutus4vars(data, countries), plutus5vars(data, countries), countries=countries, respvar="credit2gdp", shockvar="strate", filename="IR_credit2gdp_strate")
