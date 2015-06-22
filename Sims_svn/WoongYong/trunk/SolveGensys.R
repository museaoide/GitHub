SolveGensys <- function(param,SysEq,bp0,div=-1){
# returns the output from gensys for FTmodel.
# Input
#  param         : vector of parameters (should have parameter names. to be
#                  modified later
#  SysEq         : has the system of the equations of a model and necessary
#                  lists of variables (SysEq$vars, SysEq$vards, SysEq$shocks) 
#                  and the indicies of the expectational eqs (SysEq$expeqs).
#  bp0           : bond at t=0, log price level at t=0
#  div           : div for gensysct
# Output
#  fout          : gensys output with steady state values

	#---------- input for g0g1
	
	SSvars <- matrix(0,nrow=11,1)
	names(SSvars) <- SysEq$vars
	  SSvars["r"]     <- param["rbar"]
	  SSvars["a"]     <- param["abar"]
	  SSvars["p"]     <- bp0[2]								#---------- shouldn't matter
	  SSvars["w"]     <- param["pibar"]
	  SSvars["c"]     <- param["cbar"]
	  #SSvars["cc"]   <- 0
	  SSvars["b"]     <- param["taubar"]/param["thet"]
	  SSvars["tau"]   <- param["taubar"]
	  SSvars["lam"]   <- -param["sig"]*param["cbar"]
	  #SSvars["xir"]  <- 0
	  #SSvars["xipc"] <- 0
	  
	SSvards <- matrix(0,nrow=11,1)
	names(SSvards) <- SysEq$vards
		SSvards["pdot"] <- param["pibar"]
	SSshocks <- matrix(0,4,1)
	names(SSshocks) <- SysEq$shocks
	
	#---------- linearize
	
	g0g1out  <- g0g1dc(ex=SysEq$eqs,x=SysEq$vars,xd=SysEq$vards,shock=SysEq$shocks)
	g0g1oute <- g0g1evalc(dexpr=g0g1out,x=SSvars,xd=SSvards,shock=SSshocks,experr=SysEq$expeqs,param=param)
#browser()
	# normalize some eqs to avoid numerical inaccuracy
#	g0g1oute$g0["termstruc",] <- g0g1oute$g0["termstruc",]/SSvars["r"]
#	g0g1oute$g1["termstruc",] <- g0g1oute$g1["termstruc",]/SSvars["r"]
#	g0g1oute$Psi["termstruc",] <- g0g1oute$Psi["termstruc",]/SSvars["r"]
#	g0g1oute$g0["GBC",] <- g0g1oute$g0["GBC",]/SSvars["b"]
#	g0g1oute$g1["GBC",] <- g0g1oute$g1["GBC",]/SSvars["b"]
#	g0g1oute$Psi["GBC",] <- g0g1oute$Psi["GBC",]/SSvars["b"]
#	g0g1oute$g0["phcurve",] <- g0g1oute$g0["phcurve",]/SSvars["w"]
#	g0g1oute$g1["phcurve",] <- g0g1oute$g1["phcurve",]/SSvars["w"]
#	g0g1oute$Psi["phcurve",] <- g0g1oute$Psi["phcurve",]/SSvars["w"]
	g0g1oute$g0["lamdef",] <- g0g1oute$g0["lamdef",]*exp(-SSvars["lam"])
	g0g1oute$g1["lamdef",] <- g0g1oute$g1["lamdef",]*exp(-SSvars["lam"])
	g0g1oute$Psi["lamdef",] <- g0g1oute$Psi["lamdef",]*exp(-SSvars["lam"])

	# constant term
  c0 <- g0g1oute$g0 %*% SSvards - g0g1oute$g1 %*% SSvars
#browser()
  # normalize
#  c0[4] <- c0[4]/SSvars["w"]
	#dimnames(c0) <- list(names(eqs),"c")
	
	#---------- gensys

	fout <- gensysct(g0=g0g1oute$g0, g1=g0g1oute$g1, c0=c0, psi=g0g1oute$Psi, pi=g0g1oute$Pi, div=div)

  #---------- existence/uniqueness and penalty for explosive roots
  foutpnlty <- list(pnlty=0)
	if (fout$eu[2] == 1) {
	  if (fout$eu[1] == 0)	{
	  	foutpnlty <- Pnlty(g0=g0g1oute$g0, g1=g0g1oute$g1, c0=c0, psi=g0g1oute$Psi, pi=g0g1oute$Pi, div=div)
	  	# when newdiv>0, run gensys with new div
	  	if (foutpnlty$newdiv>div) {
	  	  cat("\n")
	  	  cat("adjust div and add penalty       ----------:\n")
		  	cat("# of unstable roots = ",foutpnlty$nunstab,"\n")
	  	  cat("newdiv  = ",foutpnlty$newdiv,"\n")
	  	  cat("penalty = ",foutpnlty$pnlty,"\n")
	  		fout <- gensysct(g0=g0g1oute$g0, g1=g0g1oute$g1, c0=c0, psi=g0g1oute$Psi, pi=g0g1oute$Pi, div=foutpnlty$newdiv)
	  	} else {
	  	  cat("\n")
		  	cat("non existence (far from boundary)----------: ",fout$eu,"\n")
		  	cat("# of unstable roots = ",foutpnlty$nunstab,"\n")
		  }
	  }
  }
	fout$pnlty <- foutpnlty

	fout$SSvars <- SSvars
	fout$SSvards <- SSvards
	fout$SSshocks <- SSshocks

	return(fout)
} ## THE END