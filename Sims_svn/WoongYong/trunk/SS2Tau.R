 SSvars <- matrix(0,nrow=13,1)
 names(SSvars) <- SysEq$vars
 SSvars["r"]     <- param["rbar"]
 SSvars["a"]     <- param["abar"]
 SSvars["p"]     <- bp0[2]              #---------- shouldn't matter
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
