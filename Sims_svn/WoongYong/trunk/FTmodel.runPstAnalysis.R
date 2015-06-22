i## clear workspace
rm(list=ls())
options(error = recover)
options(scipen = 5)

## user setting

machine_ 			<- 3				  # 1: Laptop on C, 2: Laptop on H, 3: Server
data_ 				<- 3					# 1: full sample 1954Q4~2005Q1, 2: 1954Q4~1979Q2, 3: 1983Q1~2005Q1
tryno_        <- 8
dev_ 					<- 2          # 1: X11(), 2: pdf()
computeCuts_ 	<- 1          # 1: compute posterior density cuts 0: don't

div_          <- 0.001      # growth rate restrictions for gensysct

# define files

logfile_ 			<- paste("./result/",data_,"/FTmodel.",data_,".xhat.",tryno_,".PstAnalysis.log",sep="")
xhatfile_			<- paste("./result/",data_,"/FTmodel.",data_,".xhat.",tryno_,".RData",sep="")
pdfname_      <- paste("./result/",data_,"/FTmodel.",data_,".xhat.",tryno_,".",sep="")

cat("=== with bond value data ===\n")
cat(sprintf("Started at     : %s \n",Sys.time()))
cat(sprintf("subsample      :  %s \n",data_))
cat(sprintf("try number     : %2.0f \n",tryno_))
cat("\n")
cat(sprintf("div            : %9.8f \n",div_))

## set working directory
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

## load sources
source("FTmodel.runLibraries.R")
source(file.path(libpath,"VARtools/trunk/rfvar3.R"))

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

nobs <- dim(y)[1]
nx <- dim(y)[2]

#---------- load parameter estimates

load(xhatfile_)
#estparam <- TrfParamBack(pstout.new$xh)
estparam <- TrfParamBack(pstout.optim$par)

estparam["thet"] <- estparam["rho"] - (1-estparam["sig"])*estparam["gamm"]
estparam["abar"] <- estparam["gamm"] + estparam["thet"] + estparam["pibar"]
estparam["rbar"] <- estparam["abar"]

#----------
pstout2 <- PstDnsty7(param=estparam,SysEq=SysEq,z=y,bp0=bp0,trf=F,dfactor=1,div=div_,prior.param=prior.param)

#---------- 
# basic output
#---------

sink(logfile_,append=FALSE,type="output")#-------------------------------

cat("\n-------------------------------------------------------------------\n")

cat("Parameter estimates\n")
print(estparam)
cat(sprintf("log prior             %18.8f\n",sum(pstout2$prd)))
cat(sprintf("log likelihood        %18.8f\n",pstout2$lkh))
cat(sprintf("log posterior density %18.8f\n",pstout2$lkh+sum(pstout2$prd)))

cat("\n-------------------------------------------------------------------\n")

sink(type="output")#----------------------------------------------------


#---------- 
# plot plots
#---------

paramexpr <- expression(phi0=phi[0],
              phi1=phi[1],
              phi2=phi[2],
              gamm=gamma,
              rho=rho,
              sig=sigma,
              bet=beta,
              delt=delta,
              omega=omega,
              psi=psi,
              pibar=bar(pi),
              cbar=bar(c),
              taubar=bar(tau),
              rhor=rho[r],
              rhopc=rho[PC],
              sig2m=sigma[M]^2,
              sig2tau=sigma[tau]^2,
              sig2r=sigma[r]^2,
              sig2pc=sigma[PC]^2,
              rhob=rho[b],
              sig2b=sigma[b]^2,
              thet=theta,
              abar=bar(a),
              rbar=bar(r))
              
obsexpr <- expression(r=r,p=p,c=c,b=b,tau=tau)
varexpr <- expression(r=r,a=a,p=p,w=w,c=c,cc=dot(c),b=b,tau=tau,lam=lambda,xir=xi[r],xipc=xi[PC],xib=xi[b])

#--------- impulse responses

fout <- SolveGensys(param=estparam,SysEq=SysEq,bp0=bp0,div=div_)

irfout <- impulsct(impact=fout$impact%*%diag(sqrt(estparam[c("sig2m","sig2tau","sig2r","sig2pc")])), G1=fout$G1, interval=1, span=40)
dimnames(irfout) <- list(SysEq$vars,SysEq$shocks,NULL)

if (dev_==2) pdf(file=paste(pdfname_,"irfs.pdf",sep="")) else X11();

		plot(ts(t(irfout[1:9, 1, ]), frequency=4, start=0), main="Responses to Monetary Shock",yax.flip=TRUE)
		plot(-ts(t(irfout[1:9, 2, ]), frequency=4, start=0), main="Responses to Fiscal Shock",yax.flip=TRUE)
		plot(ts(t(irfout[1:9, 3, ]), frequency=4, start=0), main="Responses to IS Shock",yax.flip=TRUE)
		plot(ts(t(irfout[1:9, 4, ]), frequency=4, start=0), main="Responses to PC Shock",yax.flip=TRUE)

if (dev_==2) dev.off()

#--------- one-step ahead forecasts (y(t|t-1))

yhat <- ts(matrix(0,nrow=dim(y)[1],ncol=dim(y)[2]),start=Q1,end=Q4,frequency=4)

if (dev_==2) pdf(file=paste(pdfname_,"forecastanderrors.pdf",sep="")) else X11()

		par(mfrow=c(3, 2), mar=c(2.5, 2, 2, 1), cex=.9, font=1)
		for (indx in (1:5)){
		  yhat[,indx] <- pstout2$yhat.vec[indx,,]
			plot(y[,indx],col="darkgray",main=obsexpr[indx])
			lines(yhat[,indx],col="green")
		}

		par(mfrow=c(3, 2), mar=c(2.5, 2, 2, 1), cex=.9, font=1)
		for (indx in (1:5)){
			plot(y[,indx]-yhat[,indx],main=obsexpr[indx])
		}
		#mtext(text="One-step ahead forecast errors",side=3,outer=TRUE)

if (dev_==2) dev.off()

#---------- Comparison with BVAR

start_ <- 7
vout <- rfvar3(ydata=y,lags=(start_-1),xdata=NULL,const=TRUE,breaks=NULL,lambda=5,mu=2,ic=NULL)

fcsterrors <- ts(matrix(0,nrow=nobs,ncol=nx),start=Q1,end=Q4,frequency=4)
colnames(fcsterrors) <- colnames(y)

if (dev_==2) pdf(file=paste(pdfname_,"residuals.pdf",sep="")) else X11()

		par(mfrow=c(3, 2), mar=c(4, 4, 3, 1), cex=.9, font=1)
		for (indx in (1:5)){
		  fcsterrors[,indx] <- y[,indx] - pstout2$yhat.vec[indx,,]
			plot(fcsterrors[,indx],main=paste("Residuals for",obsexpr[indx]),col="darkgray")
			lines(window(vout$u[,indx],end=Q4),col="green")
		}

		par(mfrow=c(3, 2), mar=c(4, 4, 3, 1), cex=.9, font=1)
		for (indx in (1:5)){
			sunflowerplot(x=window(fcsterrors[,indx],start=Q2,end=Q3),y=window(vout$u[,indx],start=Q2,end=Q3),main=paste("Residuals for",obsexpr[indx]),xlab="FTmodel",ylab="BVAR")
			abline(a=0,b=1,col="green")
		}

		par(mfrow=c(3, 2), mar=c(4, 4, 3, 1), cex=.9, font=1)
		for (indx in (1:5)){
			sunflowerplot(x=window(fcsterrors[,indx]^2,start=Q2,end=Q3),y=window(vout$u[,indx]^2,start=Q2,end=Q3),main=paste("Squared residuals for",obsexpr[indx]),xlab="FTmodel",ylab="BVAR")
			abline(a=0,b=1,col="green")
		}

if (dev_==2) dev.off()		

# VCV
VCV.FT <- crossprod(window(fcsterrors,start=Q2,end=Q3))/(nobs-start_+1)
VCV.BVAR <- crossprod(window(vout$u,start=Q2,end=Q3))/(nobs-start_+1)

sink(logfile_,append=TRUE,type="output")#-------------------------------

cat("\n---------------------------------------- \n");
cat("Residuals VCV: FTmodel\n");
print(VCV.FT)
cat("\n Residuals VCV: BVAR\n");
print(VCV.BVAR)
cat("\n Residuals Variances (FTmodel | BVAR)\n");
print(cbind(FTmodel=diag(VCV.FT),BVAR=diag(VCV.BVAR)))
cat("---------------------------------------- \n");

sink(type="output")#-------------------------------

#--------- latent variables

shat <- ts(matrix(0,nrow=dim(y)[1],ncol=12),start=Q1,end=Q4,frequency=4)

if (dev_==2) pdf(file=paste(pdfname_,"sshat.pdf",sep=""))  else X11()

	par(mfrow=c(4, 3), mar=c(2, 2, 1, 1), cex=.9, font=1)
	for (indx in (1:12)){
		shat[,indx] <- pstout2$shat.vec[indx+1,,]
		plot(shat[,indx],main=varexpr[indx],col="darkgray")
	}

if (dev_==2) dev.off()

#--------- likelihood cuts

if (computeCuts_==1){
	
	paramlims <- c(-.3,.3,     #phi0
	               -.3,.3,      #phi1
	               -.3,.3,      #phi2
	               -.3,.3,  #gamm
	               -.3,.3,  #rho
	               -.3,.3,        #sig
	               -.3,.3,      #bet
	               -.3,.3,#delt
	               -.3,.3,      #omega
	               -.3,.3,        #psi
	               -.3,.3,  #pibar
	               -.3,.3,        #cbar
	               -.3,.3,        #taubar
	               -.3,.3, 	#rhor
	               -.3,.3,    #rhopc
	               -.3,.3,  #sig2m
	               -.3,.3,  		#sig2tau
	               -.3,.3,  #sig2r
	               -.3,.3,  #sig2pc
	               -.3,.3,	#rhob
	               -.3,.3)	#sig2b
	paramlims <- matrix(paramlims,nrow=2,21)
	
	param1 <- estparam[1:21]
	xs <- list()
	ys <- list()
	for (indx in (1:21)){
		x1 <- seq(from=(estparam[indx]*(1+paramlims[1,indx])),to=estparam[indx],length.out=200)
		if (indx!=20) {
			x2 <- seq(from=estparam[indx],to=(estparam[indx]*(1+paramlims[2,indx])),length.out=201)
		} else {
			x2 <- seq(from=estparam[indx],to=min(.9999999999,(estparam[indx]*(1+paramlims[2,indx]))),length.out=201)
		}
		xs[[indx]] <- c(x1,x2[2:201])
		ys[[indx]] <- matrix(0,nrow=400,ncol=4)
		param2 <- param1
		for (indy in (1:400)){
			param2[indx] <- xs[[indx]][indy]
			tempout <- PstDnsty7(param=param2,SysEq=SysEq,z=y,bp0=bp0,trf=F,dfactor=1,div=div_,prior.param=prior.param)
			ys[[indx]][indy,1] <- tempout$lkh
			ys[[indx]][indy,2] <- sum(tempout$prd)
			ys[[indx]][indy,3] <- tempout$pnlty$pnlty
			ys[[indx]][indy,4] <- tempout$lkh + sum(tempout$prd)
			cat("indx: ",SysEq$params[indx]," indy: ",indy," lkh : ", ys[[indx]][indy,1]," prd: ",ys[[indx]][indy,2]," Penalty: ", ys[[indx]][indy,3]," PstDnsty: ",ys[[indx]][indy,4],"\n")
		}
	}
	
	if (dev_==2) pdf(file=paste(pdfname_,"pstcuts.pdf",sep="")) else X11()
	
			par(mfrow=c(6, 2), mar=c(2, 2, 1, 1), cex=.8)
			for (indx in (1:11)){
				plot(xs[[indx]],ys[[indx]][,4],type="l",main=paramexpr[indx])
				abline(v=estparam[indx],col="red",lwd=1.5)
			}

			par(mfrow=c(6, 2), mar=c(2, 2, 1, 1), cex=.8)
			for (indx in (12:21)){
				plot(xs[[indx]],ys[[indx]][,4],type="l",main=paramexpr[indx])
				abline(v=estparam[indx],col="red",lwd=1.5)
			}
	
	if (dev_==2) dev.off()
	
} #--------- likelihood cuts
