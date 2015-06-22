exportBISTables=function(results){
	
	#export A0 tables
	#export relLambda tables	
	
	varnames=c("Deflator", "Real GDP", "Credit-to-GDP", "Real Prop. Prices", "S.T. Rate")	
	nVars=length(varnames)
	colnames(results$A)=varnames
	rownames(results$A)=varnames
	library(xtable)
	print(xtable(results$A, caption="Rows are the response variables and the columns, the shock variables", align=rep("l", nVars+1)))

	countries=unique(results$opt$listData$countries)
	nCountries=length(countries)
	SecondaryBrk=results$opt$listData$SecondaryBrk
	nSecondaryBrks=length(SecondaryBrk)+1
	relLambda=results$relLambda

	Start=results$opt$listData$Start
	End=results$opt$listData$End

for (i in 1:nCountries){
iCountry=countries[i]
irelLambda=relLambda[, (i-1)*nSecondaryBrks + (1:nSecondaryBrks)]/relLambda[, (i-1)*nSecondaryBrks+1]

	rowNames=array()
	for (iVar in 1:nVars){
	rowNames[iVar]=paste("Shock", iVar, sep=" ") 
	}

	columnNames=array()
	for (j in 1:nSecondaryBrks){
	columnNames[j]=paste(as.yearqtr(Start[(i-1)*nSecondaryBrks + j]),as.yearqtr(End[(i-1)*nSecondaryBrks +j]), sep="-")
	}

	rownames(irelLambda)=rowNames
	colnames(irelLambda)= columnNames
	title=paste("Country ", toupper(iCountry), ":Relative variances of shocks scaled to the first period", sep="")
	
	print(xtable(irelLambda, label=toupper(iCountry), digits=rep(4, nSecondaryBrks+1), align=rep("l", nSecondaryBrks+1), caption=title))
}

}
	