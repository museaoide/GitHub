impulsdtrf <- function(B,smat,nstep)
### By:             As emerges from rfvar, neqn x nvar x lags array of rf VAR coefficients.
### smat:           nshock x nvar matrix of initial shock vectors.  To produce "orthogonalized
###                 impulse responses" it should have the property that crossprod(t(smat))=sigma,
###                 where sigma is the Var(u(t)) matrix and u(t) is the rf residual vector.  One
###                 way to get such a smat is to set smat=t(chol(sigma)).  To get the smat
###                 corresponding to a different ordering, use
###                 smat = t(chol(P %*% Sigma %*% t(P)) %*% P), where P is a permutation matrix.
###                 To get impulse responses for a structural VAR in the form A(L)y=eps, with
###                 Var(eps)=I, use B(L)=-A_0^(-1)A_+(L) (where A_+ is the coefficients on strictly
###                 positive powers of L in A), smat=A_0^(-1).
###                 In general, though, it is not required that smat be invertible.
### response:       nvar x nshocks x nstep array of impulse responses.
###
### Code written by Christopher Sims, based on 6/03 matlab code.  This version 3/27/04.
### Added dimension labeling, 8/02/04.
  {
    ##-----debug--------
    ##browser()
    ##------------------
    neq <- dim(B)[1]
    nvar <- dim(B)[2]
    lags <- dim(B)[3]
    dimnB <- dimnames(B)
    if(dim(smat)[2] != dim(B)[2]) stop("B and smat conflict on # of variables")
    response <- array(0,dim=c(neq,nvar,nstep+lags-1));
    response[ , , lags] <- smat
    response <- aperm(response, c(1,3,2))
    irhs <- 1:(lags*nvar)
    ilhs <- lags * nvar + (1:nvar)
    response <- matrix(response, ncol=neq)
    B <- B[, , seq(from=lags, to=1, by=-1)]  #reverse time index to allow matrix mult instead of loop
    B <- matrix(B,nrow=nvar)
    for (it in 1:(nstep-1)) {
      #browser()
      response[ilhs, ] <- B %*% response[irhs, ]
      irhs <- irhs + nvar
      ilhs <- ilhs + nvar
    }
    ## for (it in 2:nstep)
##       {
##         for (ilag in 1:min(lags,it-1))
##           response[,,it] <- response[,,it]+B[,,ilag] %*% response[,,it-ilag]
##       }
    dim(response) <- c(nvar, nstep + lags - 1, nvar)
    response <- aperm(response[ , -(1:(lags-1)), ], c(1, 3, 2)) #drop the zero initial conditions; array in usual format
    dimnames(response) <- list(dimnB[[1]], dimnames(smat)[[2]], NULL)
    ## dimnames(response)[2] <- dimnames(smat)[1]
    ## dimnames(response)[1] <- dimnames(B)[2]
    return(response)
  }
