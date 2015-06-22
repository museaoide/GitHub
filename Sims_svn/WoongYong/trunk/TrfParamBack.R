TrfParamBack <- function(x){

	xout <- c(x[1],                     # phi0    (-inf,inf) policy inertia
	x[2],                               # phi1    (-inf,inf) policy responsiveness on inflation
	x[3],				       									# phi2    (-inf,inf) policy responsiveness on consumption growth
	x[4],						                    # gamm    (-inf,inf) consumption growth rate
	exp(x[5]),        		              # rho	    [0,inf) time discount factor
	exp(x[6]),                          # sig     [0,inf) risk aversion
	exp(x[7]),                          # bet     [0,inf) coefficient on inflation in PC
	exp(x[8]),                          # delt    [0,inf) coefficient on consumption gap in PC
	exp(x[9])/(1+exp(x[9])),            # omega   [0,1] fiscal policy responsiveness on consumption growth
  exp(x[10]),                         # psi     [0,inf) penaly due to habit
  x[11],						                  # pibar   (-inf,inf) steady state inflation
  x[12],						                  # cbar    (-inf,inf) steady state consumption growth rate
  x[13],						                  # taubar  (-inf,inf) steady state primary surplus
  exp(x[14]),                         # rhor    [0,inf) AR coefficient on xi r (IS shock)
  exp(x[15]),                         # rhopc   [0,inf) AR coefficient on xi PC (PC shock)
  exp(x[16]),                         # sig2m   (0,inf) variance of monetary policy shock
  exp(x[17]),                         # sig2tau (0,inf) variance of fiscal policy shock
  exp(x[18]),                         # sig2r   (0,inf) variance of IS shock
  exp(x[19]),                         # sig2pc  (0,inf) variance of PC shock
  exp(x[20])/(1+exp(x[20])),					# rhob    [0,1] AR coefficient on xi b (measurement error for b)
  exp(x[21]))													# sig2b   (0,inf) variance of measurement error shock

return(xout)
}