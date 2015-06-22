## clear workspace
rm(list=ls())
options(error = recover)
options(scipen = 5)

## user setting
machine_ 			<- 3				  # 1: Laptop on C, 2: Laptop on H, 3: Server
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

estparam <- c(phi0=0.5,
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
            rhob=0.95,
            sig2b=1000)

estparam["thet"] <- estparam["rho"] - (1-estparam["sig"])*estparam["gamm"]
estparam["rbar"] <- estparam["abar"] <- estparam["gamm"] + estparam["thet"] + estparam["pibar"]

fout <- SolveGensys(param=estparam,SysEq=SysEq,bp0=bp0,div=0.001)
 
  Gkf <- matrix(0,nrow=13,ncol=13)
  Gkf[1,1]       <- 1
  Gkf[2:12,2:12] <- padm(fout$G1)
  Gkf[2:12,1]    <- numIntCnst(A=-fout$G1,m=30,delta=1)%*%fout$C
  Gkf[13,13]     <- estparam["rhob"]  #--- measurement error for b is included after solving the model

	Omg <- numIntCov(A=-fout$G1,C=fout$impact%*%diag(estparam[c("sig2m","sig2tau","sig2r","sig2pc")])%*%t(fout$impact),m=30,delta=1)
	Omg <- (Omg+t(Omg))/2
	eigOmg <- eigen(Omg)$values
	if (any(eigOmg<0)) {
		Omg <- Omg + diag(-eigOmg[11]*10,nrow=11)
	}
	MM <- matrix(0,nrow=13,ncol=13)
	MM[2:12,2:12] <- Omg
	MM[13,13] <- estparam["sig2b"]

# Cholesky decomposition for MM
M <- chol(MM[2:13,2:13])

x <- matrix(0,nrow=201,ncol=13)
colnames(x) <- c("cnst",SysEq$vars,"xib")
x[1,] <- c(cnst=1,fout$SSvars,xib=0)
for (indt in (2:201)) {	
	u <- c(0,t(M) %*% rnorm(n=12))
	x[indt,] <- Gkf %*% x[indt-1,] + u
}

x <- x[2:201,]

for (indt in (1:200)) {
	x[indt,"c"] <- x[indt,"c"] + estparam["gamm"]*indt
	x[indt,"tau"] <- x[indt,"tau"]*exp(estparam["gamm"]*indt)
  x[indt,"b"] <- x[indt,"b"]*exp(estparam["gamm"]*indt)
}

X11()
plot(ts(x[,1:6]))
X11()
plot(ts(x[,7:13]))

  Hkf <- matrix(0,nrow=5,ncol=13)
  Hkf[1,2] <- Hkf[2,4] <- Hkf[3,6] <- Hkf[4,8] <- Hkf[4,13] <- Hkf[5,9] <- 1

y <- t(Hkf %*% t(x))
colnames(y) <- c("r","p","c","b","tau")

X11()
plot(ts(y))

y2 <- y
#save(y2,estparam,file="FTmodel.fakedata.RData")