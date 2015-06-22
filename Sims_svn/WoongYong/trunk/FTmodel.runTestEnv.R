## clear workspace
rm(list=ls())
options(error = recover)
options(scipen = 5)

## user setting
machine_ 			<- 2				  # 1: Laptop on C, 2: Laptop on H, 3: Server
data_ 				<- 1					# 1: full sample 1954Q4~2005Q1, 2: 1954Q4~1979Q2, 3: 1983Q1~2005Q1
div_          <- 0.001      # growth rate restrictions for gensysct

#----------
# define files, set working directory, load sources
#----------

## load sources
if (machine_ == 1){
		setwd("C:/work/RA/Sims_2007/FTmodel3")
		libpath <- "C:/work/lib"
		dyn.load("C:/work/lib/LAPACK/blas.dll")
		dyn.load("C:/work/lib/LAPACK/lapack.dll")
} else {
	if (machine_ == 2){
			setwd("H:/work/RA/Sims_2007/FTmodel3")
			libpath <- "C:/work/lib"
			dyn.load("C:/work/lib/LAPACK/blas.dll")
			dyn.load("C:/work/lib/LAPACK/lapack.dll")
	} else {
	 libpath <- "~/work/lib"
	 dyn.load("/usr/lib/liblapack.so")
	}
}

## load libraries
source("FTmodel.runLibraries.R")

#---------- 
# load data and system of solution equations
#---------- 

load("FTmodel.data.RData")
y <- y2
load("FTmodel.SysEq.RData")

## data selection (1: full sample 1954Q4~2005Q1, 2: 1954Q4~1979Q2, 3: 1983Q1~2005Q1)
if (data_==1) {
	source("FTmodel.prior.101.R")
	bp0 <- c(830.9533,2.920847)
	} else {
	if (data_==2) {
		source("FTmodel.prior.102.R")
		y <- window(y,end=c(1979,2))
		bp0 <- c(830.9533,2.920847)
	} else {
		source("FTmodel.prior.103.R")
		y <- window(y,start=c(1983,1))
		bp0 <- c(1412.998,4.110415)
	}
}

#---------- 
# 
#---------- 

#load("./result/1/FTmodel.1.xhat.8.RData")
#load("./result.old2/2/FTmodel.2.xhat.6.RData")
#rm(pstout.new,pstout.old)

#estparam <- TrfParamBack(pstout.optim$par)
#estparam["thet"] <- estparam["rho"] - (1-estparam["sig"])*estparam["gamm"]
#estparam["rbar"] <- estparam["abar"] <- estparam["gamm"] + estparam["thet"] + estparam["pibar"]

#fout <- SolveGensys(param=estparam,SysEq=SysEq,bp0=bp0,div=0.01)
#
#irfout <- impulsct(impact=fout$impact%*%diag(sqrt(estparam[c("sig2m","sig2tau","sig2r","sig2pc")])), G1=fout$G1, interval=1, span=40)
#dimnames(irfout) <- list(SysEq$vars,SysEq$shocks,NULL)
#
#plot(ts(t(irfout[1:9, 1, ]), frequency=4, start=0), main="Responses to Monetary Shock",yax.flip=TRUE)
#plot(-ts(t(irfout[1:9, 2, ]), frequency=4, start=0), main="Responses to Fiscal Shock",yax.flip=TRUE)

#pstd <- PstDnsty(param=estparam,SysEq=SysEq,z=y,bp0=bp0,trf=FALSE,dfactor=1,div=div_,prior.param=prior.param)