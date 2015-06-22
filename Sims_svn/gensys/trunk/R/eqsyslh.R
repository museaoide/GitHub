eqsyslh <- function(eq, param, y, v0, sigvec,logged, H){
  ## eq is an equation system that has already been passed through g0g1d so that
  ##    its equation expressions have gradient expression attributes.
  ## v0 is a guess at a starting value for the steady state calculation
  ## y is the data matrix, possibly including dummy observations
  ## sigvec is the standard deviations of the structural shocks
  ## logged is a vector of names of variables that should be logged.
  ## H is the observation matrix, with dimnames.  Usually a selection matrix.
  ## It is assumed that the model is stationary.  If there are a priori
  ## certain unit roots, they should be removed from the system by appropriate differencing.
  ##------- Find ss ---------------
  ssout <- ssSolve(eq, v0, param, crit=1.07, itmax=100, verbose=FALSE)
  if(ssout$csout$rc != 0) {
    cat("No steady state.")
    if(csout$rc <4)
      cat(" Numerical problems. rc=", rc, "\n")
    else
      cat("itmax hit.\n")
    llh <- -1e20
  } else {
    vbar <- ssout$xss
    eqg0g1 <- g0g1eval(eq, x=vbar, xl = vbar)
    gout <- with(eqg0g1, gensys(g0,g1, psi=Psi, pi=Pi))
    nv <- dim(gout$G1[1])
    nsh <- dim(gout$Psi)[2]
    if (!identical(gout$eu, c(1,1))) {
      if(identical(gout$eu[1],0)) cat("Non-existence.\n")
      if(identical(gout$eu[2],0)) cat("Non-uniqueness.\n")
      llh <- -1e20
    } else {
      ## expand system to add constant, cbar
      ## logged <- c("pi","cd","cbg","ccd") # Not b, tau because these could flip sign
      vbar[logged,1] <- log(vbar[logged,1])
      ## H (observation weights) matrix.  Must be inside lh evaluation, because its elements depend on ss, and thus
      ## on param.
       ny <- dim(y)[2]
      H <- matrix(0, ny, nv+2)
      ## data are r, pi, cd + cbar, b/cd, a, tau/cd.  Note that b and tau data are normalized by c (i.e. gdp)
      ## Since cd (in original, unlogged system) is formed as c/cbar, b, the ratio of debt to cbar, divided by cd, is
      ## debt divided by c (i.e. debt over gdp in this simple model).  So data on primary surplus and on debt as
      ## ratios to gdp are used.
      xvnames <- c(dimnames(gout$G1)[[1]], "cbar", "ones")
      dimnames(H) <- list(dimnames(y)[[2]], xvnames)
      H[1, "r"] <- 1
      H[2, "pi"] <- 1
      H[3, "cd"] <- 1
      H[3, "cbar"] <- 1
      H[4, "b"] <- exp(-vbar["cd"])
      H[4, "cd"] <- -vbar["b"] * exp(-vbar["cd"])
      H[5,"a"] <- 1
      H[6, "tau"] <- exp(-vbar["cd"])
      H[6, "cd"] <- -vbar["tau"] * exp(-vbar["cd"])
      G <- cbind(gout$G1, matrix(0, nv, 2))
      G <- rbind(G, matrix(0, 2, nv+2))
      G[(nv + 1):(nv + 2), (nv + 1):(nv + 2)] <- matrix(c(1, parm["mu"],  0,1), 2, 2, byrow=TRUE)
      impact <- rbind(gout$impact, matrix(0, 2, nsh))
      impact[nv+1, "epsg"] <- 1
      ## Pull out steady states from stationary part.
      yss <- H %*% rbind(vbar, matrix(0,2,1))
      nt <- dim(y)[1]
      y <- y - matrix(yss, nt,ny, byrow=TRUE)
      ## convert G, impact, to log form
      sc <- matrix(1, nv+2, 1)
      sc[logged, ] <- vbar[logged, ]
      G <- ( (1/sc) %*% t(sc) ) * G
      impact <- (1/sc) %*% impact
      ## Prior on initial state.  Flat on cbar, SS distn for stationary part
      ssvcv <- diag(nv+2)
      ssvcv[1:nv,1:nv] <- doubling(G[1:nv,1:nv], impact[1:nv, ] %*% t(impact[1:nv, ]))
      sscv[nv+2,nv+2] <- 0
      sscv[nv+1,nv+1] <- 1e16 #This is arbitrary, to make the prior "near flat" on initial cbar. Nixes model comparisons.
      ## Now KF through the data.
      vbar <- cbind(vbar, matrix(c(0,1),2,1))
      shat <- vbar
      sig <- ssvcv
      llh <- c(0,0)
      for (it in 1:nt) {
        kout <- kf2(t(y[it, ]), H, shat, sig, G, t(impact))
        shat <- kout$shatnew
        sig <- kout$signew
        lh <- lh + kout$lh
      }
    }
    return(sum(lh))
  }
}
