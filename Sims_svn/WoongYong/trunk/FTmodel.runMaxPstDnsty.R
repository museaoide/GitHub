## clear workspace
rm(list=ls())
options(error = recover)
options(scipen = 5)
ptm <- proc.time()
sink(type="output")

## user setting
machine_ 			<- 3				  # 1: Laptop on C, 2: Laptop on H, 3: Server
data_ 				<- 3					# 1: full sample 1954Q4~2005Q1, 2: 1954Q4~1979Q2, 3: 1983Q1~2005Q1
tryno_        <- 8
output_       <- 2          # 1: output to screen, 2: output to file

div_          <- 0.001      # growth rate restriction for gensysct

#----------
# define files, set working directory, load sources
#----------

# define files

logfile_ 			<- paste("./result/",data_,"/FTmodel.",data_,".xhat.",tryno_,".log",sep="")
xhatfile_			<- paste("./result/",data_,"/FTmodel.",data_,".xhat.",tryno_,".RData",sep="")

if (file.exists(logfile_)) {
	ANSWER <- readline(prompt = sprintf("Log file %s exists. Delete and proceed? (Enter for YES)",logfile_))
	if (ANSWER=="") file.remove(logfile_) else stop(sprintf("Check log file %s",logfile_))
}
if (file.exists(xhatfile_)) {
	ANSWER <- readline(prompt = sprintf("xhat file %s exists. Delete and proceed? (Enter for YES)",xhatfile_))
	if (ANSWER=="") file.remove(xhatfile_) else stop(sprintf("Check log file %s",xhatfile_))
}

cat("=== with bond value data ===\n")
cat(sprintf("Started at     : %s \n",Sys.time()))
cat(sprintf("subsample      :  %s \n",data_))
cat(sprintf("try number     : %2.0f \n",tryno_))
cat("\n")
cat(sprintf("div            : %9.8f \n",div_))

#-------------------------------
if (output_==2) sink(logfile_,append=FALSE,type="output")
#-------------------------------

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
	 dyn.load("/usr/lib/liblapack.so.3")
	}
}

## load libraries
source("FTmodel.runLibraries.R")

#---------- 
# load data and system of solution equations
#---------- 

load("FTmodel.data.RData")
#load("FTmodel.fakedata.RData")
y <- y2 # choose data
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
# optimize
#---------- 

param0 <- c(phi0=0.5,
            phi1=0.8,
            phi2=1.5,
            gamm=0.01,
            rho=0.001,
            sig=2,
            bet=0.1,
            delt=2,
            omega=0.8,
            psi=2,
            pibar=0.01,
            cbar=7.65,
            taubar= 6,
            rhor=0.1,
            rhopc= 0.1,
            sig2m=0.000007,
            sig2tau=1.5,
            sig2r=0.000005,
            sig2pc=0.000007,
            rhob=0.9,
            sig2b=100)

#param0 <- c(phi0=.1,
#            phi1=.2,
#            phi2=1.5,
#            gamm=0.01,
#            rho=0.002,
#            sig=3,
#            bet=0.2,
#            delt= 0.2,
#            omega=0.001,
#            psi=4,
#            pibar=0.01,
#            cbar=7.65,
#            taubar= 0.06,
#            rhor=0.2,
#            rhopc= 0.2,
#            sig2m=0.00001,
#            sig2tau=0.015,
#            sig2r=0.00001,
#            sig2pc=0.000007,
#            rhob=0.9,
#            sig2b=100)

param0["cbar"] <- switch(data_,7.65,7.65,8.57)

#load("./result/2/FTmodel.2.xhat.5.RData")
#param0 <- TrfParamBack(pstout.optim$par)
#param0["phi1"] <- param0["phi1"]*.6
#param0["phi2"] <- param0["phi2"]*.6
#param0["pibar"] <- 0.0185
#param0["rho"] <- 0.0007
#param0["beta"] <- 0.03
#param0["rhopc"] <- 0.17
#param0["psi"] <- 1.7
#param0["taubar"] <- 4
#param0["sig2b"] <- 1600
#rm(pstout.new,pstout.old)

cat("\n========================================================\n")

cat(sprintf("Started at     : %s \n",Sys.time()))
cat(sprintf("subsample      :  %s \n",data_))
cat(sprintf("try number     : %2.0f \n",tryno_))
cat("\n")
cat(sprintf("div            : %6.5f \n",div_))

cat("\n--------------------------------------------------------\n")

cat("Prior used\n")
print(prior.param)

cat("\n--------------------------------------------------------\n")

cat("Initial values\n")
print(param0, digits=8)
pstd0 <- PstDnsty(param=param0,SysEq=SysEq,z=y,bp0=bp0,trf=FALSE,dfactor=1,div=div_,prior.param=prior.param)
cat("\n--------------------------------------------------------\n")
cat(sprintf("log prior             %18.8f\n",sum(pstd0$prd)))
cat(sprintf("log likelihood        %18.8f\n",pstd0$lkh))
cat(sprintf("log posterior density %18.8f\n",pstd0$lkh+sum(pstd0$prd)))
cat("--------------------------------------------------------\n")

incrmt <- 1
trf.param0 <- TrfParam(param0)
pstout.old <- list(fh=WrapPstDnsty(trf.param0,SysEq=SysEq,z=y,bp0=bp0,trf=TRUE,dfactor=-1,div=div_,prior.param=prior.param))
H0 <- diag(21)*1e-4
while (incrmt > .Machine$double.eps*10){
	pstout.new <- csminwel(fcn=WrapPstDnsty,x0=trf.param0,H0=H0,SysEq=SysEq,z=y,bp0=bp0,trf=TRUE,dfactor=-1,div=div_,prior.param=prior.param,nit=10000,crit=.Machine$double.eps*10)
	incrmt <- pstout.old$fh - pstout.new$fh
	
	cat("\n********************************************\n")
	cat("Improvement in posterior density: ", incrmt,"\n")
	print(TrfParamBack(pstout.new$xh));
	cat("\n New optimization:\n")
	cat("********************************************\n")
	
	trf.param0 <- pstout.new$xh
	pstout.old <- pstout.new
	save(pstout.old,pstout.new,param0,file=xhatfile_)
}

# once again with derivative-free algorithm
pstout.optim <- optim(par=pstout.old$xh, fn=WrapPstDnsty, gr = NULL, SysEq=SysEq,z=y,bp0=bp0,trf=TRUE,dfactor=-1,div=div_,prior.param=prior.param,method = "Nelder-Mead",control = list(trace=10,maxit=100000,reltol=.Machine$double.eps*10), hessian = FALSE)

xhat <- pstout.optim$par
save(pstout.old,pstout.new,pstout.optim,param0,file=xhatfile_)

#---------- output

param0["thet"] <- param0["rho"] - (1-param0["sig"])*param0["gamm"]
param0["abar"] <- param0["gamm"] + param0["thet"] + param0["pibar"]
param0["rbar"] <- param0["abar"]

estparam <- TrfParamBack(xhat)
estparam["thet"] <- estparam["rho"] - (1-estparam["sig"])*estparam["gamm"]
estparam["abar"] <- estparam["gamm"] + estparam["thet"] + estparam["pibar"]
estparam["rbar"] <- estparam["abar"]

cat("\n--------------------------------------------------------\n")

cat("optimization result\n")
print(estparam, digits=4)
pstd <- PstDnsty(param=estparam,SysEq=SysEq,z=y,bp0=bp0,trf=FALSE,dfactor=1,div=div_,prior.param=prior.param)
cat("\n--------------------------------------------------------\n")
cat(sprintf("log prior             %18.8f\n",sum(pstd$prd)))
cat(sprintf("log likelihood        %18.8f\n",pstd$lkh))
cat(sprintf("log posterior density %18.8f\n",pstd$lkh+sum(pstd$prd)))
cat(sprintf("penalty               %18.8f\n",pstd$pnlty$pnlty))

cat("\n--------------------------------------------------------\n")

cat(sprintf("Log file       : %s\n",logfile_))
cat(sprintf("xhat file      : %s\n",xhatfile_))
cat(sprintf("Started at     : %s \n",Sys.time()))
elapsedtime <- proc.time() - ptm
cat(sprintf("Elapsed (s)    : %15.1f\n",elapsedtime[3]))
cat(sprintf("Elapsed (h)    :           %5.1f\n",elapsedtime[3]/3600))

cat("--------------------------------------------------------\n")

#-------------------------------
if (output_==2) sink(type="output")
#-------------------------------

H <- numHess2(fcn=WrapPstDnsty,x=estparam[1:21],delta=NULL,SysEq=SysEq,z=y,bp0=bp0,trf=FALSE,dfactor=1,div=div_,prior.param=prior.param)
S <- solve(-H)

cat("diagonal elements of Hessian at mode\n")
print(diag(H))
cat("diagonal elements of inverse Hessian at mode\n")
print(diag(S))
cat("eigenvalues of inverse Hessian at mode\n")
print(eigen(S)$values)