TrfParam <- function(y){

	yout <- c(y[1],                     # phi0    (-inf,inf) policy inertia
	y[2],                               # phi1    (-inf,inf) policy responsiveness on inflation
	y[3],				       									# phi2    (-inf,inf) policy responsiveness on consumption growth
	y[4],						                    # gamm    (-inf,inf) consumption growth rate
	log(y[5]),        		              # rho	    [0,inf) time discount factor
	log(y[6]),                          # sig     [0,inf) risk aversion
	log(y[7]),                          # bet     [0,inf) coefficient on inflation in PC
	log(y[8]),                          # delt    [0,inf) coefficient on consumption gap in PC
	log(y[9]/(1-y[9])),                 # omega   [0,1] fiscal policy responsiveness on consumption growth
  log(y[10]),                         # psi     [0,inf) penaly due to habit
  y[11],						                  # pibar   (-inf,inf) steady state inflation
  y[12],						                  # cbar    (-inf,inf) steady state consumption growth rate
  y[13],						                  # taubar  (-inf,inf) steady state primary surplus
  log(y[14]),                         # rhor    [0,inf) AR coefficient on xi r (IS shock)
  log(y[15]),                         # rhopc   [0,inf) AR coefficient on xi PC (PC shock)
  log(y[16]),                         # sig2m   (0,inf) variance of monetary policy shock
  log(y[17]),                         # sig2tau (0,inf) variance of fiscal policy shock
  log(y[18]),                         # sig2r   (0,inf) variance of IS shock
  log(y[19]),                         # sig2pc  (0,inf) variance of PC shock
  log(y[20]/(1-y[20])),								# rhob    [0,1] AR coefficient on xi b (measurement error for b)
  log(y[21]))													# sig2b   (0,inf) variance of measurement error shock

return(yout)
}