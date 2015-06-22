## clear workspace
rm(list=ls())
options(error = recover)
options(scipen = 5)

## user setting

machine_ 			<- 3				  # 1: Laptop on C, 2: Laptop on H, 3: Server
data_ 				<- 3					# 1: full sample 1954Q4~2005Q1, 2: 1954Q4~1979Q2, 3: 1983Q1~2005Q1

# define files

logfile_ 			<- paste("./result/",data_,"/FTmodel.",data_,".DescStat.log",sep="")

## set working directory
if (machine_ == 1){
		setwd("C:/work/RA/Sims_2007/FTmodel2")
		libpath <- "C:/work/lib"
		dyn.load("C:/work/lib/LAPACK/blas.dll")
		dyn.load("C:/work/lib/LAPACK/lapack.dll")
} else {
	if (machine_ == 2){
			setwd("H:/work/RA/Sims_2007/FTmodel2")
			libpath <- "C:/work/lib"
			dyn.load("C:/work/lib/LAPACK/blas.dll")
			dyn.load("C:/work/lib/LAPACK/lapack.dll")
	} else {
	 libpath <- "~/work/lib"
	 dyn.load("/usr/lib/liblapack.so")
	}
}

## load sources
source("FTmodel.runLibraries.R")

#---------- load data and system of solution equations
load("FTmodel.data.RData")
y <- y2
load("FTmodel.SysEq.RData")

## data selection (1: full sample 1954Q4~2005Q1, 2: 1954Q4~1979Q2, 3: 1983Q1~2005Q1)
if (data_==1) {
	source("FTmodel.prior.101.R")
	bp0 <- c(830.9533,2.920847)
	Q1 <- c(1954,4)
	Q4 <- c(2005,1)
	Q2 <- c(1956,2)
	Q3 <- c(2005,1)
} else {
	if (data_==2) {
		source("FTmodel.prior.102.R")
		y <- window(y,end=c(1979,2))
		bp0 <- c(830.9533,2.920847)
		Q1 <- c(1954,4)
		Q4 <- c(1979,2)
		Q2 <- c(1956,2)
		Q3 <- Q4
	} else {
		source("FTmodel.prior.103.R")
		y <- window(y,start=c(1983,1))
		bp0 <- c(1412.998,4.110415)
		Q1 <- c(1983,1)
		Q4 <- c(2005,1)
		Q2 <- c(1984,3)
		Q3 <- Q4
	}
}

# descriptive statistics

capture.output(cat("Descrptive statistics \n"),file=logfile_,append=FALSE)
capture.output(summary(y),file=logfile_,append=TRUE)
capture.output(cat("\n standard deviations: \n"),file=logfile_,append=TRUE)
capture.output(sd(y),file=logfile_,append=TRUE)

dy <- y - lag(y,-1)
capture.output(cat("\n Inflation rate and consumption growth: \n"),file=logfile_,append=TRUE)
capture.output(summary(dy[,2:3]),file=logfile_,append=TRUE)
capture.output(cat("\n standard deviations: \n"),file=logfile_,append=TRUE)
capture.output(sd(dy[,2:3]),file=logfile_,append=TRUE)