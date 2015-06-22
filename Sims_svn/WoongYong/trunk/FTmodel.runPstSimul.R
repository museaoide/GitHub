## clear workspace
rm(list=ls())
options(error = recover)
options(scipen=5)

## user setting
machine_ <- 2				# 1: Laptop, 2: Server
data_ <- 2					# 1: full sample 1954Q4~2005Q1, 2: 1954Q4~1979Q2, 3: 1983Q1~2005Q1

#### setting ####

# file names
setnfile <- c(start.value="xhat4.RData",		# starting value and mean of independence MH
              save.result="xhat4.pstSimul.1.RData")	# where result is saved

# number of draws
nsimul <- 100000

# print every int.print draws
int.print <- 20

# save every int.draw draws
int.draw <- 5000

# constant
cc0 <- 1e-10
cc <- 1e-10

############################################################################

## load sources
if (machine_ == 1){
		setwd("H:/work/RA/Sims_2007/FTmodel")
		libpath <- "C:/work/lib"
		dyn.load("C:/work/lib/LAPACK/blas.dll")
		dyn.load("C:/work/lib/LAPACK/lapack.dll")
} else {
	 libpath <- "~/work/lib"
}

## load sources/libraries

source("FTmodel.runLibraries.R")
library(MASS)

#---------- load data and system of solution equations
load("FTmodel.data.RData")
load("FTmodel.SysEq.RData")

y[,4] <- y[,4]/100

#---------- load posterior mode
load(setnfile["start.value"])
rm(pstout,prior.param,param0,H)

## data selection (1: full sample 1954Q4~2005Q1, 2: 1954Q4~1979Q2, 3: 1983Q1~2005Q1)
if (data_==1) {
	source("FTmodel.prior.101.R")
	bp0 <- c(830.9533/100,2.920847)
} else {
	if (data_==2) {
		source("FTmodel.prior.102.R")
		y <- window(y,end=c(1979,2))
		bp0 <- c(830.9533/100,2.920847)
	} else {
		source("FTmodel.prior.103.R")
		y <- window(y,start=c(1983,1))
		bp0 <- c(1412.998/100,4.110415)
	}
}

#######################################################################

## time
ptm <- proc.time()

sampl.old <- estparam[SysEq$params]

eigS <- eigen(S)$values
if (any(eigS < 0))	Sig <- S - diag(19)*min(eigS)

#Sig <- diag(c(.5,1,1.5,.005,.04,1.7,.05,1,.3,1,.02,.1,.16,.8,.8,5e-7,.005,3e-7,3e-7))

sampl    <- matrix(0,nrow=nsimul,ncol=22)
  colnames(sampl) <- c(SysEq$params,"thet","abar","rbar")
pstdnsty <- matrix(0,nrow=nsimul,ncol=1)
accpt <- matrix(0,nrow=nsimul,ncol=1)

#---------- draw starting point

pst.new <- -Inf
while (pst.new == -Inf) {
	sampl.new <- mvrnorm(n = 1, mu=sampl.old, Sigma=cc0*Sig, tol = 1e-6, empirical = FALSE)
	pst.new <- WrapPstDnsty(param=sampl.new,SysEq=SysEq,z=y,bp0=bp0,trf=FALSE,dfactor=1,prior.param=prior.param)
}
		sampl.old <- sampl.new
		pst.old <- pst.new
		sampl[1,1:19] <- sampl.new
		sampl[1,20] <- sampl.new[5] - (1-sampl.new[6])*sampl.new[4]
		sampl[1,22] <- sampl[1,21] <- sampl.new[4] + sampl[1,20] + sampl.new[11]
		pstdnsty[1] <- pst.new

cat("\n-----------------------------------------------------------------\n")
cat(" starting point\n")
cat("-----------------------------------------------------------------\n")
print(as.matrix(sampl.new))
cat(sprintf("\n density at starting point = %10.5f\n",pst.new))

#---------- simulation

for (indx in 2:nsimul){

	u <- runif(1, min=0, max=1)
	
  # random walk Metropolis
	sampl.new <- mvrnorm(n = 1, mu=sampl.old, Sigma=cc*Sig, tol = 1e-6, empirical = FALSE)
	pst.new <- WrapPstDnsty(param=sampl.new,SysEq=SysEq,z=y,bp0=bp0,trf=FALSE,dfactor=1,prior.param=prior.param)
	
	r <- min(exp(pst.new-pst.old),1) # acceptance probability
	accpt[indx] <- (u < r)

	# print every -th simulation
	if (indx %% int.print == 0){
	  cat("\n-----------------------------------------------------------------\n")
		cat("simulation = ",indx,"\n")
		cat("-----------------------------------------------------------------\n")
		print(as.matrix(sampl.new),digits=6)
		cat(sprintf("\n density new= %10.5f \n",pst.new))
		cat(sprintf("\n density old= %10.5f \n",pst.old))
		cat(sprintf("\n r= %5.4f, u= %5.4f, Accept? %s (%s/%s=%4.3f)\n",r, u, accpt[indx],sum(accpt),indx,(sum(accpt)/indx)))
	}

	if (accpt[indx]) {
		sampl.old <- sampl.new
		pst.old <- pst.new
		sampl[indx,1:19] <- sampl.new
		sampl[indx,20] <- sampl.new[5] - (1-sampl.new[6])*sampl.new[4]
		sampl[indx,22] <- sampl[indx,21] <- sampl.new[4] + sampl[indx,20] + sampl.new[11]
		pstdnsty[indx] <- pst.new
	} else {
		sampl[indx,] <- sampl[indx-1,]
		pstdnsty[indx] <- pstdnsty[indx-1,]
	}

	# save intermediate results
	if (indx %% int.draw == 0){
	pstsimul <- list(indx=indx,sampl=sampl,accpt=accpt,pstdnsty=pstdnsty)
	save(pstsimul,file=setnfile["save.result"])
	rm(pstsimul)
	}
}

# save the final result
	pstsimul <- list(indx=indx,sampl=sampl,accpt=accpt,pstdnsty=pstdnsty)
	save(pstsimul,file=setnfile["save.result"])

## time
cat("\n***** time elapsed: ",proc.time() - ptm,"\n")