# Priors for the parameters

prior.param <- list()

# phi0    (-inf,inf) policy inertia    
#         See the note

prior.param$phi0 <- c(alpha=1,beta=5,location=0,scale=1)

# phi1    (-inf,inf) policy responsiveness on inflation
#         N(mu,sd)

prior.param$phi1 <- c(mu=1,sd=1) ######

# phi2    (-inf,inf) policy responsiveness on consumption growth  
#         N(mu,sd)

prior.param$phi2 <- c(mu=.5,sd=1.5)

# gamm    (-inf,inf) consumption growth rate                      
#         N(mu,sd)

prior.param$gamm <- c(mu=.01,sd=.005)

# rho	    [0,inf) time discount factor                            
#         Gamma(shape,rate)

prior.param$rho <- c(shape=1.5,rate=30)

# sig     [0,inf) risk aversion
#         Gamma(shape,rate)

#prior.param$sig <- c(shape=3,rate=1)
prior.param$sig <- c(shape=4,rate=1)
                           
# bet     [0,inf) coefficient on inflation in PC
#         Gamma(shape,rate)

prior.param$bet <- c(shape=3,rate=20)

# delt    [0,inf) coefficient on consumption gap in PC
#         Gamma(shape,rate)

prior.param$delt <- c(shape=4,rate=2)

# omega   [0,1] fiscal policy responsiveness on consumption growth
#         Beta(shape1,shape2)
## adjusted

prior.param$omega <- c(shape1=1,shape2=19)

# psi     [0,inf) penaly due to habit
#         Gamma(shape,rate)

prior.param$psi <- c(shape=1,rate=1)
                        
# pibar   (-inf,inf) steady state inflation                       
#         N(mu,sd)

prior.param$pibar <- c(mu=.01,sd=.005)
#prior.param$pibar <- c(mu=.005,sd=.001)

# cbar    (-inf,inf) steady state consumption growth rate
#         N(mu,sd)
#         log(real C) = 8.55 at 1982Q4

prior.param$cbar <- c(mu=8.55,sd=.1)

# taubar  (-inf,inf) steady state primary surplus        
#         N(mu,sd)
#         price level at 1954Q3 =  exp(2.920847)
#         nominal Mkt value of debt at 1954Q3 = 154.2
#         real Mkt value of debt at 1954Q3 = 154.2/exp(2.920847)*100 = 830.9533
## adjusted

prior.param$taubar <- c(mu=.02*830.9533,sd=.02*830.9533)
         
# rhor    [0,inf) AR coefficient on xi r (IS shock)
#         Gamma(shape,rate)

prior.param$rhor <- c(shape=1.5,rate=1.5)

# rhopc   [0,inf) AR coefficient on xi PC (PC shock)
#         Gamma(shape,rate)

prior.param$rhopc <- c(shape=1.5,rate=1.5)

# sig2m   (0,inf) variance of monetary policy shock
#         IG(shape,rate)
#         rate in IG = scale in Gamma

prior.param$sig2m <- c(shape=2,rate=1e+6)

# sig2tau (0,inf) variance of fiscal policy shock                 
#         IG(shape,rate)
#         rate in IG = scale in Gamma
## adjusted

#prior.param$sig2tau <- c(shape=2,rate=.03)
prior.param$sig2tau <- c(shape=2,rate=1)

# sig2r   (0,inf) variance of IS shock
#         IG(shape,rate)
#         rate in IG = scale in Gamma

prior.param$sig2r <- c(shape=2,rate=1e+6)
                        
# sig2pc  (0,inf) variance of PC shock
#         IG(shape,rate)
#         rate in IG = scale in Gamma

prior.param$sig2pc <- c(shape=2,rate=1e+6)

# rhob    [0,1] autocorelation coefficient for measurement error
#         Beta(shape1,shape2)

prior.param$rhob <- c(shape1=8,shape2=2)

# sig2b   (0,inf) variance of measurement error shock
#         IG(shape,rate)
#         rate in IG = scale in Gamma

prior.param$sig2b <- c(shape=2,rate=.005)