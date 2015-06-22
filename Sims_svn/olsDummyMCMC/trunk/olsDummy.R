olsDummy <- function(y=NULL,X=NULL,yd=NULL,Xd=NULL,lambda=1,ndraw=0)
  ## y:             dependent variables data                                               
  ## X:             rhs variable data matrix                                              
  ## yd:            dependent variables for dummy observations scaled by sigma             
  ## Xd:            rhs variable matrix with yd
  ## lambda:        scale factor for dummy obs
  ## Note:          There is no prior specified here on Sigma=Var(y|X), as this can be done through
  ##                dummy observations.  If y is Txm, then p dummy observations with the Xd matrix
  ##                0 and the yd matrix r generates an inverse-Wishart prior on Sigma with
  ##                scale matrix r'r and degrees of freedom p-m-1.  If you want a different degrees of freedom
  ##                parameter, increase (or decrease) T by the amount of increased or decreased degrees of freedom
  ##                you want in the Sigma prior.  [The same dummy observation, interpreted as generating a prior on
  ##                the precision Sigma^(-1) rather than on Sigma, implies degrees of freedom p+m+1.  But this program
  ##                integrates with respect to Sigma, not its inverse, so the p-m-1 degrees of freedom is implied.  Also,
  ##                the program by default includes a Jeffreys-type prior factor of det(Sigma)^(-m-1)/2, so that the p
  ##                dummy observations together with the Jeffreys component imply a prior with p degrees of freedom on Sigma^(-1).
  ##                In the common case where 1) the X matrix is Txk, 2) there are q dummy observations characterizing the
  ##                distribution of the regression coefficients, 3) these q dummy observations satisfy yd=Xd %*% bd for some vector bd
  ##                (i.e., they do not make conflicting claims about the location of bd), and 4) there are p dummy observations with
  ##                yd2=r, Xd2=0, the prior implies a marginal pdf for Sigma^(-1) with scale matrix r'r and degrees of freedom p+q-k.
  ##                A proper prior therefore requires k+m dummy observations in total, generally.
  ## wfull:         log of integral of posterior pdf
  ## wp:            log of integral of prior pdf
  ## fit:           return value of the fit (via lm.fit()) of the model with actual and dummy data
  ## pfit:          return value of the fit (via lm.fit()) of the model with dummy data alone.  Allows construction of prior.
  ## cSig:          cholesky decomposition of the Sigma matrix for ndraw draws from the marginal posterior.
  ## b:             m x k x ndraw array of draws of the parameter matrix (or vector, when m=1)
  {
    m <- dim(y)[2]
    k <- dim(X)[2]
    Xs <- rbind(X,lambda*Xd)
    ys <- rbind(y,lambda*yd)
    Ts <- dim(ys)[1]
    Td <- dim(Xd)[1]
    Xd <- lambda*Xd
    yd <- lambda*yd
    ## fit <- lm.fit(y=ys,x=Xs)
    fit <- lm(ys~Xs+0)
    ## pfit <- lm.fit(y=yd,x=Xd)
    if(Td >= m+k)
      {
        pfit <- lm(yd~Xd+0)
        csp <- qr.R(qr(pfit$residuals))
        csp <- matrix(sign(diag(csp)),m,m)*csp #(flip signs of any rows with negative diagonals)
        cxp <- qr.R(pfit$qr)
        cxp <- matrix(sign(diag(cxp)),k,k)*cxp
        cxp <- t(solve(cxp))
        wp <- matrictint(T=dim(Xd)[1],cs=csp,cx=cxp)
      }
    else
      {
        wp <- NULL
      }
    cs <- qr.R(qr(fit$residuals))
    cs <- matrix(sign(diag(cs)),m,m)*cs 
    cx <- qr.R(fit$qr)
    cx <- matrix(sign(diag(cx)),k,k)*cx #(flip signs of any rows with negative diagonals)
    cx <- t(solve(cx))
    wfull <- matrictint(T=Ts,cs=cs,cx=cx)
    if(ndraw>0)
      {
        csi <- t(solve(cs))
        cSig <- rvecwish(v=Ts-k,cs=csi,ndraw=ndraw) # m x m x ndraw array or 1 x ndraw vector when m=1
        b <- rnorm(ndraw*k*m)
        dim(b) <- c(k,m*ndraw)
        b <- crossprod(cx,b)
        if (m==1)
          {
            cSig <- 1.0/cSig
            b <- cSig*b
          }
        else
          {
            for (id in 1:ndraw){cSig[,,ndraw] <- t(solve(cSig))}
            dim(b) <- c(k,m,ndraw)
            b <- aperm(b,c(2,1,3))
            for (id in 1:ndraw) {b[,,id] <- crossprod(cSig[,,id],b[,,id])}
            aperm(b,c(2,1,3))
            dimnames(b) <- list(eqn=character(1:m),vble=character(1:k),draw=NULL)
          }
      }
    else
      {
        b <- NULL
        cSig <- NULL
      }
    b <- b+array(fit$coefficients,dim=dim(b))
    return(list(w=wfull-wp,wp=wp,fit=fit,pfit=pfit,cSig=cSig,b=b))
  }        
