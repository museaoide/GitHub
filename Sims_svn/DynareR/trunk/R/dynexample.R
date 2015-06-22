source("dynread.r")

## example of using these functions
model <- c("(c/c(1))^gamma*beta*(alpha*exp(a(1))*k^(alpha-1)+1-delta)=1",
"a=rho*a(-1)+eps",
"k+c=exp(a)*k(-1)^alpha+(1-delta)*k(-1)")

dynfw <- makedynframework(dynpath="/home/tpapp/bin/dynare++",
                          var=c("k","c","a"),
                          varexo="eps",
                          parameters=c("beta", "gamma", "rho",
                            "alpha", "delta"),
                          model=model,
                          initval=c(k=0.066,c=0.43,a=0.01),
                          order=2,
                          statevars=c("a","k","eps"))

dd <- calldynare(dynfw, c(.99,          # beta
                    2,                  # gamma
                    .9,                 # rho
                    .3,                 # alpha
                    .025                # delta
                    ), matrix(0.001,1,1))

unfoldtensor(dd$dyn$dyn.g.2[1,],k=2,n=3) # 2nd order tensor for the
                                        # first variable

## match variables names
pp <- matchvarnames(dynfw$statevars,dd$dyn$dyn.state.vars)
dd$dyn$dyn.state.vars[pp]
