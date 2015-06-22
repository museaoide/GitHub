PstDnsty <- function(param,SysEq,z,bp0,trf,dfactor=1,div,prior.param){
# Input
#  param         : vector of parameters (should have parameter names. to be
#                  modified later)
#  SysEq         : has the system of the equations of a model and necessary
#                  lists of variables (SysEq$vars, SysEq$vards, SysEq$shocks) 
#                  and the indicies of the expectational eqs (SysEq$expeqs).
#  z             : data [r p c b tau]
#  bp0           : bond at t=0, log price level at t=0 (with bond value data, bp0[1] no need)
#  trf           : parameters are transformed for optimization?
#  dfactor       : factor to be multiplied to log likelihood and log prior
#                  density. e.g. -1 for minimization.
#  prior.param   : the parameters of prior distributions. to get the log likelihood
#                  only, set prior=NULL.
# Output
#

  nz <- dim(z)[1]

  #---------- when the parameters are reparameterized
  if (trf) {
  	param <- TrfParamBack(param)
  #---------- when not
  } else {
	  if (any(param[5]<0,param[6]<0,param[7]<0,param[8]<0,param[9]<0,param[9]>1,param[10]<0,param[14]<0,param[15]<0,param[16]<0,param[17]<0,param[18]<0,param[19]<0,param[20]<0,param[20]>1,param[21]<0)) {
	    cat("parameter out of the domain \n")
	    return(fv=list(lkh=0,prd=-dfactor*Inf,eu=c(-9,-9)))
	  }
 }

  # a is actually log(a) and abar = gamma + theta + pibar.
  param["thet"] <- param["rho"] - (1-param["sig"])*param["gamm"]
  param["rbar"] <- param["abar"] <- param["gamm"] + param["thet"] + param["pibar"]

	# for the consumer optimization to be well-difined (theta>0)
	if (param["thet"]<=0) {
	  cat("theta < 0:\n")
		print(param[c("rho","sig","gamm","thet")])
		return(list(lkh=0,prd=-dfactor*Inf,eu=c(-8,-8)))
	}

  #---------- 
  # gensys
  #---------- 

  fout <- SolveGensys(param=param,SysEq=SysEq,bp0=bp0,div=div)

  #---------- existence/uniqueness and penalty for explosive roots
	if (fout$eu[2] == 0) {
	  cat("\n")
  	cat("non uniqueness                   ----------: ",fout$eu,"\n")
  	return(fv=list(lkh=0,prd=-dfactor*Inf,eu=fout$eu,fout=fout,pnlty=fout$pnlty))
  } else {
	  if (fout$eu[1] == 0) {
	    if (fout$pnlty$newdiv==0)	{
	  	  # when newdiv=0, no penalty. just abandon.
		    return(fv=list(lkh=0,prd=-dfactor*Inf,eu=fout$eu,fout=fout,pnlty=fout$pnlty))
	    }
	  }
  }

  if (fout$eu[1] == -2 || fout$eu[2] == -2) return(fv=list(lkh=0,prd=-dfactor*Inf,eu=fout$eu,fout=fout,pnlty=fout$pnlty))

  #----------
  # Kalman filter
  #---------- 
        # z(t+1) = [r, p, c, b, tau]
        # s(t+1) = [1,r,a,p,w,c,cc,b,tau,lam,xir,xipc,xib]
        # z(t+1) = Hs(t+1)
        # s(t+1) = Gs(t) + Me(t+1)

  Hkf <- matrix(0,nrow=5,ncol=13)
  Hkf[1,2] <- Hkf[2,4] <- Hkf[3,6] <- Hkf[4,8] <- Hkf[4,13] <- Hkf[5,9] <- 1
  
  Gkf <- matrix(0,nrow=13,ncol=13)
  Gkf[1,1]       <- 1
  Gkf[2:12,2:12] <- padm(fout$G1)
  Gkf[2:12,1]    <- numIntCnst(A=-fout$G1,m=30,delta=1)%*%fout$C
  Gkf[13,13]     <- param["rhob"]  #--- measurement error for b is included after solving the model

  #Mkf <- matrix(0,nrow=12,ncol=4)
  #Omg <- dtInnovCov(A=-fout$G1,C=fout$impact%*%diag(param[c("sig2m","sig2tau","sig2r","sig2pc")])%*%t(fout$impact))
  #Mkf[2:12,] <- t(chol(Omg$X))
  #Mkft <- t(Mkf)
	Omg <- numIntCov(A=-fout$G1,C=fout$impact%*%diag(param[c("sig2m","sig2tau","sig2r","sig2pc")])%*%t(fout$impact),m=30,delta=1)
	Omg <- (Omg+t(Omg))/2
	eigOmg <- eigen(Omg)$values
	if (any(eigOmg<0)) {
		Omg <- Omg + diag(-eigOmg[11]*10,nrow=11)
	}
	MM <- matrix(0,nrow=13,ncol=13)
	MM[2:12,2:12] <- Omg
	MM[13,13] <- param["sig2b"]

  #---------- initial conditions
  #           1) ---
  #           2) log(price level) N(2.920847,0.1^2)
  #           3) real Mkt value of debt N(830.9533,10^2)

#  shat.old <- c(cnst=1,fout$SSvars,xib=0)
#  sig.old <- diag(c(0,r=.25*1e-4,a=.25*1e-4,p=.01,w=1e-4,c=1e-2,cc=1e-4,b=100,tau=10,lam=1,xir=.5*param["sig2r"]/param["rhor"],xipc=.5*param["sig2pc"]/param["rhopc"],xib=param["sig2b"]/(1-param["rhob"]^2)))
#browser()
  A <- matrix(0,nrow=12,ncol=12)
  A[2:12,2:12] <- fout$G1
  A[2:12,1] <- fout$C
  Omega <- matrix(0,nrow=12,ncol=12)
  Omega[2:12,2:12] <- fout$impact%*%diag(param[c("sig2m","sig2tau","sig2r","sig2pc")])%*%t(fout$impact)
	mu0 <- c(cnst=1,fout$SSvars[c("p","c","b","tau","lam")])
  ssndx <- c(1,4,6,8,9,10)
  Sig0 <- diag(1e0,nrow=12)
  Sig0[1,1] <- 0
  s0 <- SigInit(A=A, Omega=Omega, T=nz, mu0=mu0, Sig0=Sig0,ct=TRUE, ssndx=ssndx)
	shat.old <- c(Re(s0$mu),0)
	sig.old <- matrix(0,nrow=13,ncol=13)
	sig.old[1:12,1:12] <- Re(s0$v)
	sig.old[13,13] <- param["sig2b"]/(1-param["rhob"]^2)

  #---------- run kf

#browser()
  lkh <- as.numeric(-param["gamm"]*nz*(nz+1)/2)              # Jacobian for tau
  for (indx in 1:nz){
		z[indx,"c"] <- z[indx,"c"] - param["gamm"]*indx         # adjust for growth
		z[indx,"tau"] <- z[indx,"tau"]*exp(-param["gamm"]*indx)
		z[indx,"b"] <- z[indx,"b"]*exp(-param["gamm"]*indx)
    #kfout <- kf2(y=z[indx,],H=Hkf,shat=shat.old,sig=sig.old,G=Gkf,MM=MM)
    kfout <- kf(y=z[indx,],H=Hkf,shat=shat.old,sig=sig.old,G=Gkf,MM=MM)
    shat.old <- kfout$shatnew;
    sig.old <- kfout$signew;
    lkh <- lkh + kfout$lh[1] + kfout$lh[2]
  }
  lkh <- lkh - fout$pnlty$pnlty
#browser()
  #---------- 
  # prior densities
  #---------- 
	
  if (is.null(prior.param)) prd <- NULL else prd <- PriorDnsty(param,prior.param)

	if (is.na(lkh)) lkh <- -Inf
	if (any(is.na(prd))) prd <- -Inf

  return(list(lkh=dfactor*lkh,prd=dfactor*prd,eu=fout$eu,pnlty=fout$pnlty))
} ## THE END