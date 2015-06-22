sysmat <- function(B) {
## Constructs the lags*nv x lags*nv system matrix for the "stacked" first-order
## version, from the nv x nv x lags array of coefficients By returned by rfvar3.
    n <- dim(B)
    dim(B) <- c(n[1],n[2]*n[3])
    B <- rbind(B,diag(1,nrow=n[1]*(n[3]-1),ncol=n[2]*n[3]))
    return(B)
  }



trimts <- function(y) {
  if(!is.ts(y)) error("arg to trimts() not a time series")
  if(is.null(dim(y))) dim(y) <- c(length(y), 1)
  firstGood <- match(TRUE, apply(y, MARGIN=1, function(x) !any(is.na(x))))
  firstGood <- time(y)[firstGood]
  y <- window(y, start=firstGood)
  lastGood <- match(TRUE, apply(y, MARGIN=1, function(x) any(is.na(x))))
  if(!is.na(lastGood)) {
    lastGood <- time(y) [lastGood-1]
    y <- window(y, end=lastGood)
  }
  return(y)
}

rfvar3 <- function(ydata=NA,lags=nLags,xdata=NULL,nCountries=nCountries, const=const,breaks=NULL,lambda=5,mu=2,ic=NULL, sigpar=NULL) {
    ## This algorithm goes for accuracy without worrying about memory requirements.
    ## ---------------------------------------------------------------------------
    ## The standard prior it implements is NOT APPROPRIATE for seasonally unadjusted data, even
    ## if seasonal dummies are included in xdata.  The prior shrinks toward simple persistence, so it
    ## will tend to prevent the dummies from picking up all the seasonality.
    ## ---------------------------------------------------------------------------
    ## ydata:   T x nvar dependent variable data matrix.  
    ## xdata:   T x nx exogenous variable data matrix.  
    ##          Note that if either ydata or xdata has only one column, it must still have a dim vector.  In
    ##          other words it must be a Tx1 array, not a vector of length T.
    ##------------------
    ## const:   If TRUE, a column of ones is added to (or becomes, if xdata is NULL) the xdata matrix.
    ## lags:    number of lags
    ## breaks:  rows in ydata and xdata after which there is a break.  This allows for
    ##          discontinuities in the data (e.g. war years) and for the possibility of
    ##          adding dummy observations to implement a prior.  This must be a column vector.
    ##          Note that a single dummy observation becomes lags+1 rows of the data matrix,
    ##          with a break separating it from the rest of the data.  The function treats the 
    ##          first lags observations at the top and after each "break" in ydata and xdata as
    ##          initial conditions.  If not null, breaks vector[1] must equal 0.
    ## lambda:  weight on "co-persistence" prior dummy observations.  This expresses
    ##          belief that when all variables are at a fixed initial level, they tend to
    ##          stay there.  This is consistent with stationarity and with nonstationarity with or
    ##          without cointegration.  With lambda < 0 , the 
    ##          constant term is not included in the dummy observation, so that stationary models
    ##          with means equal to initial ybar do not fit the prior mean.  With lambda>0, the prior
    ##          implies that large constants are unlikely if unit roots are present.  To omit this type of
    ##          dummy observation, use lambda=NULL.
    ## mu:      weight on "own persistence" prior dummy observation.  Expresses belief
    ##          that when y_i has been stable at its initial level, it will tend to persist
    ##          at that level, regardless of the values of other variables.  There is
    ##          one of these for each variable.  A reasonable first guess is mu=2.
    ##          To omit this type of dummy observation, use mu=NULL
    ## ic:      for direct input of the initial conditions mean that is used in the persistence dummy observations,
    ##          as ic$ybar and ic$xbar. 
    ##          If is.null(ic), the mean of the first lags observations in ydata, xdata are used.
    ##      The program assumes that the first lags rows of ydata and xdata are real data, not dummies.
    ##      Dummy observations should go at the end, if any.  If pre-sample x's are not available,
    ##      repeating the initial xdata(lags+1,:) row or copying xdata(lags+1:2*lags,:) into 
    ##      xdata(1:lags,:) are reasonable subsititutes.  These values are used in forming the
    ##      persistence priors.
    ## sigpar: list(A0, lmd, Tsigbrk) Allow SVAR with time varying shock variances. Tsigbrk, if not null, must begin with 0. See below.
    ## returns:
    ## By:      nvar x nvar x lags matrix of coefficients on lagged y's.  1st dimension is "equation number"
    ## Bx:      nvar x nx matrix of coefficients on x's
    ## u:       (T-6+ (number of dummy obs)) x nvar matrix of residuals.  If ydata is a ts object, u will be also, and will
    ##          be correctly dated.  u observations dated after end(ydata) are dummy observations.
    ## xxi:     X'X inverse, same for all equations.  kronecker(cov(u),xxi) is the full covariance matrix of the regression coefficients.
    ## snglty:  Usually 0.  If the rhs variable matrix is not full column rank, this is the gap between the number of columns and the
    ##          number of non-zero singular values.
    ## Code written by Christopher Sims.  This version 8/13/04.
    ## 12/18/05:  added ts properties for u, better comments.
    ##---------------------------------------
    ## Modified 2013.8.12 to allow use of A0, lmd, Tsigbrk.  With non-null A0, By is A+ from
    ## A0 %*% y(t) = A+(L) %*% y(t) + exp(.5 lmd(t)) * eps(t) .  This works even with
    ## lmd constant, but in that case running a single rf estimate (A0=I), then iterating
    ## on (A0, lmd) alone makes more sense. With lmd varying, rf estimates change with lmd.
    ## --------------------------------------------------------
    if (is.null(dim(ydata))) dim(ydata) <- c(length(ydata),1)
    T <-dim(ydata)[1]
    ## Note that if rfvar3() has been called with dummy obs's already in place, this T
    ## includes the dummies.
    nvar<-dim(ydata)[2]
    ##nox=isempty(xdata)
    if (const) {
        xdata <- cbind(xdata,matrix(1,T,1))
    }
    nox <- identical(xdata,NULL)
    if(!nox){
        T2 <- dim(xdata)[1]
        nx <- dim(xdata)[2]
    } else {
        T2 <- T; nx <- 0; xdata<- matrix(0,T2,0)
    } 
    ## note that x must be same length as y, even though first part of x will not be used.
    ## This is so that the lags parameter can be changed without reshaping the xdata matrix.
    ## ------------------------
    if (!identical(T2,T)) {
        print('Mismatch of x and y data lengths')
        return()
    }
    if (identical(breaks,NULL))
        nbreaks <- 0
    else {
        ## if (is.ts(ydata)) {                # Can use Yr, month-or-quarter pairs, or real number dates.
        ##   if (is.matrix(breaks) ) {
        ##     breaks <- breaks[ , 1] + (breaks[ ,2] - 1) / frequency(ydata)
        ##   } else {
        ##     if (any(abs(breaks - round(breaks))) > 1e-8) {
        ##       breaks <- match(breaks, time(ydata))
        ##     }
        ##   }                               #if not real numbers, not yr-month pairs, it's just obs number
        ## }
        ## Any use of tsp(ydata) has to be in external processing functions.
        nbreaks<-length(breaks)
    }
    
    breaks <- c(0, breaks,T)
    if(any(breaks[2:length(breaks)] < breaks[1:(length(breaks)-1)]))
        stop("list of breaks must be in increasing order\n")
    ## initialize smpl as null if initial observations are only there for lambda/mu prior.
    ## matlab code uses the fact that in matlab a:b is null if b<a, which is not true for R.
    ## if(breaks[2]>lags)
    ##   smpl <- (lags+1):breaks[2]
    ## else
    ##   smpl <- NULL
    ## if(nbreaks>0){
    ##   for (nb in 2:(nbreaks+1))
    ##     smpl <- c(smpl,(breaks[nb]+lags+1):breaks[nb+1])
    ## }
    smpl <- NULL
    for (nb in 2:(nbreaks + 2)) {
        if ( breaks[nb] > breaks[nb-1] + lags )
            smpl <- c(smpl, (breaks[nb-1] + lags + 1):breaks[nb])
    }
    ## With logic above, one can use an mts-type ydata and omit sections of it by including sequences of breaks separated by
    ## less than lags+1.  E.g. with lags=6, monthly data, breaks=rbind(c(1979,8), c(1980,2), c(1980,8), c(1980,12)) omits
    ## Sep 1979 through Dec 1981, plus 6 months after that, which are initial conditions for the next sample segment.
    Tsmpl <- length(smpl)
    X <- array(0,dim=c(Tsmpl,nvar,lags))
    for(ix in seq(along=smpl))
    #grab the lags (time, var, lag)
        X[ix,,] <- t(ydata[smpl[ix]-(1:lags),,drop=FALSE])    
        
    dim(X) <- c(Tsmpl,nvar*lags)
    X <- cbind(X, xdata[smpl,,drop=FALSE])
    y <- ydata[smpl,,drop=FALSE]
    ## Everything now set up with input data for y=Xb+e 
    ## ------------------Form persistence dummies-------------------
    if (! (is.null(lambda) & is.null(mu) ) ) {
        if(is.null(ic)) {
            ybar <- apply(as.array(ydata[1:lags,,drop=FALSE]),2,mean)
            dim(ybar) <- c(1,dim(ydata)[2])
            if (!nox) {
                xbar <- apply(array(xdata[1:lags,,drop=FALSE],dim=c(lags,dim(xdata)[2])),2,mean)
                dim(xbar)=c(1,dim(xdata)[2])
            } else {
                xbar <- NULL
            }
        } else {
            ybar <- ic$ybar
            xbar <- ic$xbar
        }
        if (!is.null(lambda)){
            if (lambda<0){
                lambda <- -lambda
                xbar <- array(0,c(1,dim(xdata)[2]))
            }
            xdum <- lambda * cbind(array(rep(ybar,lags),dim=c(1,lags*length(ybar))), xbar)
            ydum <- array(0,c(1,nvar))
            ydum[1,] <- lambda*ybar
            y <- rbind(y,ydum)
            X <- rbind(X,xdum)
        }
        if (!is.null(mu)) {
            xdum <- cbind(
                array(rep(diag(as.vector(ybar),nrow=length(ybar)),lags),
                      dim=c(dim(ybar)[2],dim(ybar)[2]*lags)),
                array(0,dim=c(nvar,dim(xdata)[2])))*mu
            ydum <- mu*diag(as.vector(ybar),nrow=length(ybar))
            X <- rbind(X,xdum)
            y <- rbind(y,ydum)
        }
    }
    if (!is.null(sigpar)) {
        Tsigbrk <- sigpar$Tsigbrk
        lmd <- sigpar$lmd
        A0 <- sigpar$A0
        if (!is.null(Tsigbrk)) {
            ## Tsigbrk <- invtime(Tsigbrk, ydata) #so Tsigbrk given as dates
            nsig <- length(Tsigbrk)
        } else {
            nsig <- 1
        }
        Tsigbrk <- c(Tsigbrk, T)
        
        lmdndx <- rep(1:nsig, times=diff(Tsigbrk))
        lmdseries <- lmd[ , lmdndx]
        if ( Tsmpl < dim(y)[1] ) {      #dummy obs formed in rfvar3
            ## Should not be combining this branch with dummy obs's from varprior()
            ## already included in ydata.
            lmdp <- apply(lmdseries[ ,smpl], 1, mean)
            lmdseries <- cbind(lmdseries[ , smpl], matrix(lmdp, nvar, dim(y)[1] - Tsmpl))
        } else {
            lmdseries <- lmdseries[ , smpl]
        }
        ## i.e., use mean of lmdseries for dummy observation weights.  Note that
        ## since lmd is logged, this is geometric mean, maybe not best. 
        nX <- dim(X)[2]
        ya0 <- y %*% t(A0)
        B <- matrix(0,  nX, nvar)
        u <- matrix(0, Tsmpl, nvar)
        uraw <- u
        xxi <- array(0, c(nX, nX, nvar))
        logdetxxi <- vector("numeric", nvar)
        snglty <- vector("numeric", nvar)
        for (iq in 1:nvar) {
   
            
   #note that the * isn't matrix multiplication
   			Xq<-X
            for (iq in 1:nvar){
            	wt <- exp(.5 * lmdseries[iq, ])
            Xq[,iq*(1:lags)] <-  wt * X[,iq*(1:lags)]
            for (iC in 1:nCountries){
            	Xq[, nvar*lags+ (iC-1)*nvar*lags+ (iq-1)*lags + (1:lags)] <- wt * X[, nvar*lags+ (iC-1)*nvar*lags+ (iq-1)*lags + (1:lags)]
            }
            }
            
            wt <- exp(.5 * lmdseries[iq, ])
            yq <- wt * ya0[ , iq]
            lso <- lsfit(Xq, yq, intercept=FALSE)
            ## intercept already in X. resids should be unit vce.
            B[ , iq] <- lso$coefficients
            ##method for solving for lsf
            Rq <- qr.R(lso$qr)
            xxi[ , , iq] <- solve(crossprod(Rq))
            u[ , iq] <- lso$residuals
            logdetxxi[iq] <- -2 * sum(log(abs(diag(Rq))))
            snglty[iq] <- (logdetxxi[iq] == -Inf)
        }
    } else {
        ## Instead of svd below, could invoke lsfit.  Faster?
        vldvr <- svd(X)
        dfx <- sum(vldvr$d > 100*.Machine$double.eps)
        di <- 1./vldvr$d[1:dfx]
        vldvr$u <- vldvr$u[, 1:dfx]
        vldvr$v <- vldvr$v[, 1:dfx]
        snglty <- dim(X)[2] - dfx
        logdetxxi <- 2 * sum(log(abs(di)))
        ##B <- vldvr$v %*% diag(di,nrow=length(di)) %*% t(vldvr$u) %*% y (line below is just more efficient)
        B <- vldvr$v %*% (di * (t(vldvr$u) %*% y))
        u <-  y-X %*% B
        xxi <-  di * t(vldvr$v)
        xxi <-  crossprod(xxi)
        uraw <- NULL       # so it won't be missing in the list of outputs
    }
    if (!is.null(tsp(ydata))) u <- ts(u, start=start(ydata)+c(0,lags),freq=frequency(ydata))
    ## dates at end of sample are for dummy obs, meaningless.  If there are other
    ## nontrivial breaks, the dates for u are also meaningless.
    ## dim(B) <-  c(nvar*lags+nx,nvar) # rhs variables, equations (this was redundant)
    By <-  B[1:(nvar*lags),]
    dim(By) <-  c(nvar,lags,nvar)       # variables, lags, equations
    By <-  aperm(By,c(3,1,2)) #equations, variables, lags to match impulsdt.m
    ## label all the output, if the data matrices had labels
    if(!is.null(dimnames(ydata)[2]))
        {
            ynames <- dimnames(ydata)[[2]]
        }else
            {
                ynames <- rep("",times=nvar)
            }
    if(!nox)
        {
            if(!is.null(dimnames(xdata)[[2]]))
                {
                    xnames <- dimnames(xdata)[[2]]
                } else {
                    xnames <- rep(" ",times=nx)
                }
        }
    dimnames(By) <- list(ynames,ynames,as.character(1:lags))
    xxinames <- c(paste(rep(ynames,lags),rep(1:lags, each=length(ynames)),sep=""),xnames)
    dimnames(xxi) <- list(xxinames,xxinames)
    if (nox)
        Bx <-  NULL
    else
        {
            Bx <-  matrix(B[nvar*lags+(1:nx),],dim(B)[2],nx)
            dimnames(Bx) <- list(ynames,xnames)
        }
### logintlh <-  matrictint(u'*u,xxi,size(X,1)-nvar-1)-.5*nvar*(nvar+1)*log(2*pi);
### Might want to create a version without the dimnames if using this in a program.
    ##------------
    ## returns some things that are not available with sigpar=NULL.  Either split to
    ## separate programs, or create alternate return lists.
    if(!is.null(sigpar)) {
       return(list(By=By, Bx=Bx, u=u, uraw=uraw, xxi= xxi, snglty=snglty, logdetxxi=logdetxxi,
                lmdseries=lmdseries, call=match.call())) #var.logintlh <-  logintlh
    } else {
      return(list(By=By, Bx=Bx, u=u, xxi= xxi, snglty=snglty, logdetxxi=logdetxxi, call=match.call()))
    }
}

restrictVAR <- function(vout, type=c("3", "KF","SVhtskd"), rmat=NULL, yzrone=NULL, xzrone=NULL,
                        const=NULL, cyzr=NULL, cxzr=NULL) {
    ## restrictions can be specified as rows of rmat, with coefficients applied to elements of By and Bx
    ## stacked as they are in xxi (and then repeated across the equation index), or they can be specified
    ## in yzrone, xzrone.  Each zero element of yzrone or xzrone generates a restriction that sets the corresponding
    ## coefficient in By or Bx to zero (or to a constant, if !is.null(const)).  Both kinds of restrictions
    ## can be non-trivial in the same call.
    ##-------------------------------------------
    ## type:     vout as from rfvar3 ("3") or as from rfvarKF ("KF")
    ## const:    the right hand side of rmat %*% coeff = const, not the constant in the var.
    ## cyzr, cxzr:  If using yzrone, xzrone with non-trivial constants, leave const=NULL and specify
    ##           constants with cyzr and cxzr
    ##------------------------------------------
    ## sc:       The Schwarz criterion rejects the restriction if the chisq value plus the sc value
    ##           is positive.  This version of the sc is scale-sensitive.  Variables with higher
    ##           variance are penalized more strongly, as with a prior that expects higher-variance
    ##           variables to explain more variance.
    ##
    ##------------------------------------------------
    ## Note 2013-3-4:  Try eliminating scale effects by converting X'X to correlation matrix
    if(type == "SVhtskd") {
        require("tensor")
        bvw <- vout
        vout <- bvw$vout$var
    }
    ncf <- dim(vout$By)[2] * dim(vout$By)[3] + dim(vout$Bx)[2]
    neq <- dim(vout$By)[1]
    ny <- dim(vout$By)[2]
    lags <- dim(vout$By)[3]
    nx <- dim(vout$Bx)[2]
    if (is.null(rmat)) {
        rmat <- matrix(0, 0, ncf *neq)
    }
    if (!is.null(yzrone)) {
        byz <- which(yzrone == 0, arr.ind=TRUE)
        nrstr <- dim(byz)[1]
        if (is.null( cyzr)) cyzr <- array(0, dim(yzrone))
        for (ir in 1:nrstr ) {
            newrow <- rep(0, neq * ncf)
            newrow[(byz[ir,1] - 1) * ncf + (byz[ir, 3] -1) * ny + byz[ir, 2]] <- 1
            rmat <- rbind(rmat,newrow)
        }
        const <- c(const, cyzr[byz])
    }
    if (!is.null(xzrone)) {
        bxz <- which(xzrone == 0, arr.ind=TRUE )
        nrstr <- dim(bxz)[1]
        if (is.null(cxzr)) cxzr <- matrix(0, neq, nx)
        for (ir in 1:nrstr)  {
            newrow <- rep(0,ncf * neq)
            newrow[(bxz[ir,1] - 1) * ncf + ny * lags + bxz[ir, 2]] <- 1
            rmat <- rbind(rmat, newrow)
        }
        const <- c(const, cxzr[bxz])
    }
    svdr <- svd(rmat)
    if (max(abs(svdr$d)) > 1e10 * min(abs(svdr$d))){
        error("restrictions not full rank")
    }
    ## Note that t(rv) spans the same space as rmat, so the restrictiosn are crossprod(v,coeffs)=gamma
    ## rv <- svdr$v    #2013.5.9
    T <- if (type == "3" || type == "SVhtskd") dim(vout$u)[1] else dim(vout$fcsterr)[1]
    if (type == "3") {
        sig <- cov(vout$u)
        svdsig <- svd(sig)
        singsig <- (max(svdsig$d) > 1e10 * min(svdsig$d))
        svdxxi <- svd(vout$xxi)
        singxxi <- (max(svdxxi$d) > 1e10 * min(svdxxi$d))
        singv <- singsig || singxxi
        if(!singv) {
            ## schwarz <- rmat %*% kronecker(svdsig$u %*% diag(1/sqrt(svdsig$d)), svdxxi$u %*% diag(1/sqrt(svdxxi$d)))
            ##schwarz <- kronecker((1/sqrt(svdsig$d)) * t(svdsig$u), (1/sqrt(svdxxi$d)) * t(svdxxi$u)) %*% rv  #2013.5.9
            ## sqrtVb <- kronecker(sqrt(svdsig$d) * t(svdsig$u), 1/sqrt(svdxxi$d)) * t(svdxxi$u)
            ## line above seems to be a mistake, since xxi is already x'x-inverse
            sqrtVb <- kronecker(sqrt(svdsig$d) * t(svdsig$u), sqrt(svdxxi$d) * t(svdxxi$u))
            dgVb <- apply(sqrtVb^2, 2, sum)
            rmatC <- rmat %*% diag(sqrt(T * dgVb))
            sqrtVbC <- sqrtVb %*% diag(1/sqrt(T * dgVb))
            lndetVb <- sum(log(svdsig$d)) * dim(vout$xxi)[1] + sum(log(svdxxi$d)) * dim(sig)[1]
            lndetVbC <- lndetVb - sum(log(dgVb * T))
        }
    } else if (type == "KF") {             #type=="KF"
        svdVb <- svd(vout$Vb)
        sqrtVb <- sqrt(diag(svdVb$d)) %*% t(svdVb$u)
        dgVb <- diag(vout$Vb)
        rmatC <- rmat %*% diag(sqrt(T * dgVb))
        sqrtVbC <- sqrtVb %*% diag(1/sqrt(T * dgVb))
        lndetVb <- sum(log(svdVb$d))
        lndetVbC <- lndetVb - sum(log(dgVb * T))
        ## schwarz <- rmat %*% svdVb$u %*% diag(1/sqrt(svdVb$d)) #below is more efficient version for large Vb
        ## schwarz <- (1/sqrt(svdVb$d)) * (t(svdVb$u) %*% rv)
    } else {                            #type="SVhtskd"
        nv <- dim(vout$By)[1]
        nX <- dim(vout$xxi)[1]
        Vb <- matrix(0, nX * nv, nX * nv)
        for (iq in 1:nv) {
            Vb[nX * (iq-1) + 1:nX, nX * (iq-1) + 1:nX] <- vout$xxi[ , , iq]
        }
        A0i <- solve(bvw$A)
        Vb <- kronecker(A0i, diag(nX)) %*% Vb %*% kronecker(t(A0i), diag(nX))
        svdVb <- svd(Vb)
        sqrtVb <- sqrt(diag(svdVb$d)) %*% t(svdVb$u)
        dgVb <- diag(Vb)
        rmatC <- rmat %*% diag(sqrt(T * dgVb))
        sqrtVbC <- sqrtVb %*% diag(1/sqrt(T *dgVb))
        lndetVb <- sum(log(svdVb$d))
        lndetVbC <- lndetVb - sum(log(dgVb * T))
    }
    svdvr <- svd(sqrtVb %*% t(rmat))
    svdvrC <- svd(sqrtVbC %*% t(rmatC)) #result == line above?
    vdim1 <- dim(svdvr$u)[1]
    svdvrp <- svd(diag(vdim1) - svdvr$u %*% t(svdvr$u), nu=vdim1 - dim(rmat)[1])
    svdvrpC <- svd(diag(vdim1) - svdvrC$u %*% t(svdvrC$u), nu=vdim1 - dim(rmat)[1])
    svdvrpuv <- svd(crossprod(svdvrp$u, t(sqrtVb))) 
    svdvrpuvC <- svd(crossprod(svdvrpC$u, t(sqrtVbC))) 
    lndetUR <- sum(log(svdvrpuv$d))
    lndetURC <- sum(log(svdvrpuvC$d))
    df <- dim(rmat)[1]
    ## schwarz <- -2 * sum(log(diag(chol(crossprod(schwarz)))))   +  df * log(2 * pi)
    schwarz <- lndetVb - 2 * lndetUR + df * log(2 * pi)
    schwarzC <- lndetVbC - 2 * lndetURC + df * log(2 * pi)
    if(is.null(const)) const <- rep(0, dim(rmat)[1])
    if(type == "SVhtskd") {
        vout$By <- tensor(A0i, vout$By, 2, 1)
        vout$Bx <- A0i %*% vout$Bx
    }
    stackedcf <- c(t(cbind(matrix(vout$By, nrow=neq), vout$Bx)))
    gap <- rmat %*% stackedcf - const
    ##svdv <- svd(rmat %*% vout$Vb %*% t(rmat))
    chstat <- (1/svdvr$d) * t(svdvr$v) %*%  gap
    chstat <- crossprod(chstat)
    return(list(chiSquared=chstat, df=df, sc=schwarz, pval=pchisq(chstat,df), sc2 = schwarz - (ncf*neq-df)*log(1 - df/(neq*ncf)), scC=schwarzC ))
}

impulsdtrf <- function(vout=NULL, smat=NULL, nstep=40, order=NULL)
### vout:           output structure from rfvar3.
###                 To use this with output from postdraw, create a dummy vout with vout$By=pout$By[ , , ,id] and provide smat=pout$smat[ , ,id]
### smat:           if !is.null(vout) and order and smat are NULL, the impulse responses will be for a cholesky decomp with variables
###                 ordered as in the input to rfvar3.  More generally, can be any set of initial
###                 values for the shocks.  To keep the shocks orthogonal and scaled properly,
###                 smat should be such that smat %*% t(smat) == crossprod(vout$u)/dim(u)[1].
###                 However, the routine works with singular smat or with smat's column dimension
###                 less than its row dimension.
### order:          To get a cholesky decomp with a different ordering, set order to an integer
###                 vector giving the desired ordering.  
### response:       nvar x nshocks x nstep array of impulse responses.
###
###                 with vout from rfvarKF, smat argument is required, since there is no vout$u.
###
### Code written by Christopher Sims, based on 6/03 matlab code.  This version 3/27/04.
### Added dimension labeling, 8/02/04.  Allow non-square smat, integrate with rfvar3 output, 4.7.10.
  {
    ##-----debug--------
    ##browser()
    ##------------------
    B <- vout$By
    neq <- dim(B)[1]
    nvar <- dim(B)[2]
    lags <- dim(B)[3]
    dimnB <- dimnames(B)
    if (is.null(smat)) {
      if (is.null(order) ) {
        order <- 1:neq
      }
      smat <- t(pchol(crossprod(vout$u)/dim(vout$u)[1], order)) # makes first shock affect all variables
    }
    nshock <- dim(smat)[2]
    if(dim(smat)[1] != dim(B)[1]) stop("B and smat conflict on # of equations") #
    response <- array(0,dim=c(neq,nshock,nstep+lags-1));
    response[ , , lags] <- smat
    response <- aperm(response, c(1,3,2))
    irhs <- 1:(lags*nvar)
    ilhs <- lags * nvar + (1:nvar)
    response <- matrix(response, ncol=nshock)
    B <- B[, , seq(from=lags, to=1, by=-1)] #reverse time index to allow matrix mult instead of loop
    B <- matrix(B,nrow=nvar)
    for (it in 1:(nstep-1)) {
      response[ilhs, ] <- B %*% response[irhs, ]
      irhs <- irhs + nvar
      ilhs <- ilhs + nvar
    }
    ## for (it in 2:nstep)
    ##       {
    ##         for (ilag in 1:min(lags,it-1))
    ##           response[,,it] <- response[,,it]+B[,,ilag] %*% response[,,it-ilag]
    ##       }
    dim(response) <- c(nvar, nstep + lags - 1, nshock)
    response <- aperm(response[ , -(1:(lags-1)), ,drop=FALSE], c(1, 3, 2)) #drop the zero initial conditions; array in usual format
    dimnames(response) <- list(dimnB[[1]], dimnames(smat)[[2]], NULL)
    ## dimnames(response)[2] <- dimnames(smat)[1]
    ## dimnames(response)[1] <- dimnames(B)[2]
    return(response)
  }


pchol <- function(sig, porder) {
  invporder <- match(1:length(porder), porder)
  return(chol(sig[porder, porder])[invporder, invporder])
}


rwwish <- function (v, S,n) {
  ## v:  degrees of freedom for the Wishart
  ## S:  Scale matrix for the Wishart
  ## n:  number of draws needed
  ## value:  a pxpxn array of n draws that are Cholesky square roots of Wisharts.
  ## The routine could easily be modified (as shown below) to deliver Wisharts
  ## instead of their square roots, but for simulation purposes usually the square
  ## roots are more useful, and the user can add the line of code to convert after
  ## the call in any case.
  ##
  ## This is much more efficient than the original MCMCpack rwish when multiple
  ## draws are to be made.  about 20 times faster for 4x4 S with n=200, e.g.
  ## Because this is based on
  ## MCMCpack code, it is covered by the GPL, version 2, which is available with
  ## every R installation.
  if (!is.matrix(S)) 
    S <- matrix(S)
  p <- nrow(S)
  if (  p != ncol(S)) {
    stop(message = "S not square.\n")
  }
  if (v < p) {
    stop(message = "df v is less than the dimension of S.\n")
  }
  CC <- chol(S)
  Z <- matrix(0, n, p^2)
  dseq <- (0:(p-1))*p+(1:p)
  nch <- n*p
  Z[ , dseq] <- sqrt(rchisq(nch, v:(v - p + 1)))
  if (p > 1) {
    pseq <- 1:(p - 1)
    Z[ , rep(p * pseq, pseq) + unlist(lapply(pseq, seq))] <-
      rnorm(n * p * (p - 1)/2)
  }
  dim(Z) <- c(n*p,p)
  Z <- Z %*% CC
  dim(Z) <- c(n, p,p)
  Z <- aperm(Z,c(2,3,1))
  ## for (id in 1:n){
  ##   Z[ , , id] <- crossprod(Z[ , , id])
  ## }                                     
  return(Z)
}


postdraw <- function(vout,n,nosigprior=FALSE){
## Provides n draws from the joint posterior on residual covariance and coefficients of a
## reduced form VAR.
## vout:	return value from rfvar3()
## n:		Number of draws
## nosigprior:	If true, don't use the Jeffreys-motivated improper prior.
##--------------------
  xxi <- chol(vout$xxi)
  ncf <- dim(xxi)[1]
  S <- crossprod(vout$u)
  df <- dim(vout$u)[1]-dim(xxi)[1]      # This effectively uses a |Sigma|^{-(nvar+1)/2} prior
  neq <- dim(vout$u)[2]
  lags <- (ncf-dim(vout$Bx)[2])/neq
  if(nosigprior){df <- df-dim(S)[1]-1}	# This undoes the |Sigma|^{-(nvar+1)/2} prior
  wmat <- rwwish(df,solve(S),n)
  for (it in 1:n){wmat[,,it] <- solve(wmat[,,it])}
  nmat <- array(rnorm(n*neq*ncf),c(ncf,neq,n))
  cfmat <- t(cbind(matrix(vout$By,neq,neq*lags),vout$Bx))
  for (ir in 1:n){
    nmat[,,ir] <- crossprod(xxi,nmat[,,ir]) %*% wmat[,,ir]+cfmat
  }
  Byx <- aperm(nmat, c(2,1,3))
  By <- Byx[ , 1:(neq*lags), ]
  dim(By) <- c(neq,neq,lags,n)
  ## Bx <- as.vector(vout$Bx)+aperm(nmat,c(2,1,3))[,(neq*lags+1):ncf,]  # Bx added in both here and in cfmat. Bug caught by A.Zawadwoski
  Bx <- Byx[,(neq*lags+1):ncf,]
  ## Note that wmat[,,i] is ready for use as input to impulsdtrf
  return(list(By=By,Bx=Bx,smat=wmat))
}


varprior <-  function(nv=1,nx=0,lags=1,mnprior=list(tight=5,decay=.5),vprior=list(sig=1,w=1),
                      urprior=list(lambda=NULL, mu=NULL), xsig=NULL, ybar=NULL, xbar=1, nstat=rep(TRUE,nv))
### ydum, xdum:   dummy observation data that implement the prior
### breaks:       vector of points in the dummy data after which new dummy obs start
###                   Set breaks=T+matrix(c(0,breaks),ncol=1), ydata=rbind(ydata,ydum), xdum=rbind(xdata,xdum), where 
###                   actual data matrix has T rows, in preparing input for rfvar3
### nv,nx,lags: VAR dimensions
### mnprior$tight:Overall tightness of Minnesota prior. 1/tight ~ own lag std dev
### mnprior$decay:Standard deviations of lags shrink as lag^(-decay)
### vprior$sig:   Vector of prior modes for square roots of diagonal elements of r.f. covariance matrix
###                  Names of this vector name columns of output ydum.
### vprior$w:     Weight on prior on vcv.  1 corresponds to "one dummy observation" weight
###                   vprior$sig is needed
###                   to scale the Minnesota prior, even if the prior on sigma is not used itself.
###                   Set vprior$w=0 to achieve this.
###                   mnprior and vprior.w can each be set to NULL, thereby eliminating the corresponding
###                   dummy observations.
### xsig:          rough scale of x variances.  names of this vector name output xdum
### urprior:       Parameters of the "unit roots" and "co-persistence" priors that are
###                   implemented directly in rfvar3.  lambda and mu should be NULL here if
###                   the dummy observations generated here are used with rfvar3 and lanbda and mu
###                   are not NULL in rfvar3.   lambda < 0 means x'st not included.  Note that constant
###                   is assumed to be last element of x.  If you want lambda < 0 to be the only source
###                   of a prior on the constant, but xsig is not null, set the last element of xsig
###                   to zero.  
### ybar,xbar:        estimates of data means, used in constructing urprior component, but not otherwise.
###                   The default xbar=1 is correct when the constant is the only x.    
### nstat:         Set components corresponding to non-persistent variables to FALSE.
### Note:          The original Minnesota prior treats own lags asymmetrically, and therefore
###                   cannot be implemented entirely with simple dummy observations.  It is also usually
###                   taken to include the sum-of-coefficients and co-persistence components
###                   that are implemented directly in rfvar3.R.  The diagonal prior on v, combined
###                   with sum-of-coefficients and co-persistence components and with the unit own-first-lag
###                   prior mean generates larger prior variances for own than for cross-effects even in 
###                   this formulation, but here there is no way to shrink toward a set of unconstrained 
###                   univariate ARs.
###-----------------------
###
{ require(abind)
  ## nx=0 case messes up, at least at the end (2012.9.23)
  if (!is.null(mnprior))
    { ## single-coefficient prior dummy obs.
      ## each vbl and each lag has a dummy observation, and each dummy obs has values for current and lagged
      ## y's  and current x's. we separate the y's and the x's into two arrays.  The last two indexes, lag
      ## and rhsy, index the dummy observations.  
      xdum <- if(nx > 0) {
        array(0, dim=c(lags + 1, nx, lags, nv), dimnames=list(obsno=1:(lags + 1), xvbl=1:nx, lag=1:lags, rhsy=1:nv))
      } else {
        NULL
      }
      ydum <- array(0,dim=c(lags+1,nv,lags,nv),dimnames=list(obsno=1:(lags+1),rhsy=1:nv,lag=1:lags, rhsy=1:nv))
      for (il in 1:lags)
        {
          ##-----debug---------
          ## browser()
          ##------------------
          ydum[il+1,,il,] <- il^mnprior$decay*diag(vprior$sig,nv,nv)
        }
      ## If we have non-trivial x's, need dobs's for them, also.
      if(!is.null(xsig)) {
        ydumx <-  array(0, dim=c(lags + 1, nv, nx), dimnames=list(obsno=1:(lags + 1), rhsy=1:nv, dx=1:nx))
        xdumx <-  array(0, dim=c(lags + 1, nx, nx), dimnames=list(obsno=1:(lags + 1), xvbl=nx, dx=1:nx))
        xdumx[1, , ] <- diag(xsig, nx, nx)
        ## note that xvalues for obsno 2:(lags+1) don't matter.  This is one dummy obseervation,
        ## so only the "current" x is used.
      }
      ydum[1,,1,] <- diag(vprior$sig * nstat, nv, nv) # so own lag has mean zero if nstat FALSE
      ydum <- mnprior$tight * ydum
      dim(ydum) <- c(lags+1,nv,lags*nv)
      ydum <- ydum[seq(lags+1,1,by=-1),,]
      xdum <- mnprior$tight*xdum
      dim(xdum) <- c(lags+1,nx,lags*nv)
      xdum <- xdum[seq(lags+1,1,by=-1),,]
    } else {
      ydum <- NULL;
      xdum <- NULL;
      breaks <- NULL;
      lbreak <- 0;
    }
  if (!is.null(urprior$lambda) ) {
    ## lambda obs.  just one
    ydumur <- matrix(ybar, nrow=lags+1, ncol=nv, byrow=TRUE) * abs(urprior$lambda)
    if(urprior$lambda > 0) {
      xdumur <- matrix(xbar, lags + 1, nx, byrow=TRUE) * urprior$lambda # (all but first row redundant)
    } else {
      xdumur <- matrix(0, lags + 1, nx)
    }
  } else {
    ydumur <- NULL
    xdumur <- NULL
  }
  ## mu obs. sum(nstat) of them
  if (!is.null(urprior$mu)) {
    ydumuri <-array(0, c(lags+1, nv, nv))
    for (iv in which(nstat)) {
      ydumuri[ , iv, iv] <- ybar[iv]
    }
    ydumur <- abind(ydumur, urprior$mu *ydumuri, along=3)
    xdumur <- abind(xdumur, array(0, c(lags+1, nx, nv)), along=3)
  }
  if (!is.null(vprior) && vprior$w > 0)
    {
      ydum2 <- array(0,dim=c(lags+1,nv,nv))
      xdum2 <- array(0,dim=c(lags+1,nx,nv))
      ydum2[lags+1,,] <- diag(vprior$sig,nv,nv)*vprior$w #The vprior$w factor was missing until 11/29/06
                                        # Original idea, not implemented, was probably that w be an integer
                                        # repetition count for variance dobs.
                                        # Now it's just a scale factor for sig. in variance prior.
    } else {
      ydum2 <- NULL
      xdum2 <- NULL
    }
  ## stack everything up.
  dim(ydum) <- c(lags + 1, nv, lags * nv) # merge all the individual mn dobs
  dim(xdum) <- c(lags + 1, nx, lags * nv)
  ydum <- abind(ydum, ydumur, ydum2, along=3)
  xdum <- abind(xdum, xdumur, xdum2, along=3)
  breaks <- (lags+1) * (1:(dim(ydum)[3] -1)) # end of sample is not a "break".
  ydum <- aperm(ydum, c(1, 3, 2))
  ydum <- matrix(ydum, ncol=dim(ydum)[3])
  xdum <- aperm(xdum, c(1,3,2))
  xdum <- matrix(xdum, ncol=dim(xdum)[3])
  ##   dim(ydum2) <- c((lags+1)*nv,nv)
  ##   dim(ydum) <- c((lags+1)*nv,lags*nv)
  ##   ydum <- cbind(ydum,ydum2)
  ##   dim(xdum2) <- c((lags+1)*nx,nv)
  ##   dim(xdum) <- c((lags +1)*nx,lags*nv)
  ##   xdum <- cbind(xdum,xdum2)
  ##   dim(ydum) <- c(lags+1,nv,dim(ydum)[2])
  ##   ydum <- aperm(ydum,c(1,3,2))
  ##   dim(ydum) <- c(dim(ydum)[1]*dim(ydum)[2],nv)
  ##   dim(xdum) <- c(lags+1,nx,dim(xdum)[2])
  ##   xdum <- aperm(xdum,c(1,3,2))
  ##   dim(xdum) <- c(dim(xdum)[1]*dim(xdum)[2],nx)
  ##   if(nv>1){
  ##     breaks <- c(breaks, (lags+1)*(0:(nv-1))+lbreak)
  ##   }
  ## } else {
  ##   if (!is.null(ydum)) { # case with mnprior non-null, but vprior null
  ##     ydum <- aperm(ydum, c(1, 3, 2))
  ##     dim(ydum) <- c(prod(dim(ydum)[1:2]), dim(ydum)[3])
  ##     xdum <- aperm(xdum, c(1,3,2))
  ##     dim(xdum) <- c(prod(dim(xdum)[1:2]), dim(xdum)[3])
  ##   }
  ## }
  dimnames(ydum) <- list(NULL, names(vprior$sig))
  dimnames(xdum) <- list(NULL, names(xsig))
  return(list(ydum=ydum,xdum=xdum,pbreaks=breaks))
  ## data here in the form of T by nv y, and T x nx x.  Lagged y's not put in to a rhs
  ## regression matrix, so a "breaks" vector is needed.  
  ## rfvar3 adds persistence and sum of coeffs dummy observations at end of  data in lhs and rhs
  ## regression matrix form.  So to combine this with rfvar3, set lambda and mu to NULL in one or the
  ## other program.
}

bfgsi <- function(H0,dg,dx) {
### dg is previous change in gradient; dx is previous change in x;
### 6/8/93 version that updates inverse hessian instead of hessian
### itself.
### Copyright by Christopher Sims 1996.  This material may be freely
### reproduced and modified.
  n <- length(dg)
  dim(dg) <- c(n,1)
  dim(dx) <- c(n,1)
  Hdg <- H0 %*% dg
  dgdx <- as.numeric(crossprod(dg,dx))
  dxHdg <- drop(dx %o% Hdg) # drops are needed to get rid of redundant unit-dimensions
  x1 <- as.numeric(1+crossprod(dg,Hdg)/dgdx)*drop(dx %o% dx)
  x2 <- dxHdg+t(dxHdg)
  ## code below causes problems when x1 and x2 have matching zeros, and I can't now (2005-12-22)
  ## figure out why I thought it was a good idea
##   if ( max(abs(x1-x2)/(abs(x1)+abs(x2))) <= 1e-12 ) {
##     cat("bfgs update failed.\n")
##     cat("x1, x2 = ",x1,x2,"\n")
##     H=H0
##     return(H)
##   }
  if (abs(dgdx) <= 1e-12)   {
    cat("bfgs update failed.\n")
    cat("|dg| =", sqrt(sum(dg^2)), "|dx| = ", sqrt(sum(dx^2)),"\n")
    cat("crossprod(dg,dx) =", dgdx,"\n")
    cat("|H*dg| =", crossprod(Hdg),"\n")
    H=H0
    return(H)
  }
  ## otherwise, 
  H <- H0 + x1/dgdx - x2/dgdx
  save(file="H.dat", H)
  return(H)
}


csminit <- function(fcn,x0,f0,g0,badg,H0,...){
### retcodes: 0, normal step.  5, largest step still improves too fast.
### 4,2 back and forth adjustment of stepsize didn't finish.  3, smallest
### stepsize still improves too slow.  6, no improvement found.  1, zero
### gradient.
###---------------------
### Fixed 7/17/93 to use inverse-hessian instead of hessian itself in bfgs
### update.
###
### Fixed 7/19/93 to flip eigenvalues of H to get better performance when
### it's not psd.
###
  ## ANGLE <- .0005
  ANGLE <- 1e-7
  THETA <- .01 #(0<THETA<.5) THETA near .5 makes long line searches, possibly fewer iterations.
  FCHANGE <- 1000
  MINLAMB <- 1e-9
### fixed 7/15/94
### MINDX <- .0001;
### MINDX <- 1e-6;
  MINDFAC <- .01
  fcount<-0
  lambda<-1
  xhat<-x0
  f<-f0
  fhat<-f0
  g <- g0
  gnorm <- sqrt(sum(g^2))
###
  if (!badg && (gnorm < 1.e-12)) {      # put !badg 8/4/94
    retcode <- 1
    dxnorm <- 0
    ## gradient convergence
  } else {
    ## with badg true, we don't try to match rate of improvement to directional
    ## derivative.  We're satisfied just to get some improvement in f.
    ##
    ##if(badg)
    ##   dx = -g*FCHANGE/(gnorm*gnorm);
    ##  dxnorm = norm(dx)
    ##  if dxnorm > 1e12
    ##     disp('Bad, small gradient problem.')
    ##     dx = dx*FCHANGE/dxnorm;
    ##   end
    ##else
    ## Gauss-Newton step;
    ##---------- Start of 7/19/93 mod ---------------
    ##[v d] = eig(H0);
    ##toc
    ##d=max(1e-10,abs(diag(d)));
    ##d=abs(diag(d));
    ##dx = -(v.*(ones(size(v,1),1)*d'))*(v'*g);
    dx <- -H0 %*% g
    dxnorm <- sqrt(sum(dx^2))
    if (dxnorm > 1e12){
      cat("Near-singular H problem.\n")
      dx <- dx*FCHANGE/dxnorm
    }
    dfhat <- crossprod(dx,g0)
    if(!badg){
      ## test for alignment of dx with gradient and fix if necessary
      a <- -dfhat/(gnorm*dxnorm)
      if(a<ANGLE){
        if (a < 0) {
          dx <- -g / gnorm^2
          dfhat <- -1
          cat("H unused\n")
          ## Don't bother with H.  It's not psd. Step length here appropriate for log LH,
          ## where 1.0 is a reasonable scale for changes.
        } else {
          dx <- dx - as.numeric(ANGLE*dxnorm/gnorm+dfhat/(gnorm*gnorm))*g
          dx <- dx * dxnorm / sqrt(sum(dx^2)) # This line added 2/18/2004
          dfhat <- crossprod(dx,g)
          ## dxnorm <- sqrt(sum(dx^2)) # No longer needed with 2/18/2004 change
          cat("Correct for low angle:" ,a,"\n")
        }
      }
    }
    cat(sprintf("Predicted improvement            = %18.9f\n",-dfhat/2))
    ## cat("Predicted improvement:", sprintf("%18.9f",-dfhat/2),"\n")
    ##
    ## Have OK dx, now adjust length of step (lambda) until min and
    ## max improvement rate criteria are met.
    done <- 0
    factor <- 3
    shrink <- 1
    lambdaMin <- 0
    lambdaMax <- Inf
    lambdaPeak <- 0
    fPeak <- f0
    lambdahat <- 0
    while(!done){
      ## argument of fcn retains its dim (or lack thereof), but g
      ## always emerges as n x 1
      ddx <- dx*lambda
      dim(ddx) <- dim(x0)
      dxtest <- x0+ddx
      f <- fcn(dxtest,...)
      cat(sprintf("lambda = %10.5g; f = %20.7f",lambda,f ),"\n")
      if(f < fhat){
        fhat <- f
        xhat <- dxtest
        lambdahat <- lambda
      }
      fcount <- fcount+1
      shrinkSignal <- (!badg && (f0-f < max(-THETA*dfhat*lambda,0))) || (badg && ((f0-f) < 0)) 
      growSignal <- (!badg && ( (lambda > 0)  &&  (f0-f > -(1-THETA)*dfhat*lambda) ))
      if(  shrinkSignal  &&   ( (lambda>lambdaPeak) || (lambda<0) )){
        if ((lambda>0) && ((!shrink || (lambda/factor <= lambdaPeak)))){
          shrink <- 1
          factor <- factor^.6
          while(lambda/factor <= lambdaPeak){
            factor <- factor^.6
          }
          if( abs(factor-1)<MINDFAC){
            if( abs(lambda) < 4){
              retcode <- 2
            }else{
              retcode <- 7
            }
            done <- 1
          }
        }
        if(  (lambda<lambdaMax) && (lambda>lambdaPeak) ){
          lambdaMax <- lambda
        }
        lambda <- lambda/factor
        if( abs(lambda) < MINLAMB ){
          if( (lambda > 0) & (f0 <= fhat) ){
            ## try going against gradient, which may be inaccurate
            lambda <- -lambda*factor^6
          }else{
            if( lambda < 0 ){
              retcode <- 6
            }else{
              retcode <- 3
            }
            done <- 1
          }
        }
      }else{
        if(  (growSignal && lambda>0) ||  (shrinkSignal && ((lambda <= lambdaPeak) && (lambda>0)))  ) {
          if( shrink ){
            shrink <- 0
            factor <- factor^.6
            if( abs(factor-1)<MINDFAC ) {
              if( abs(lambda)<4 ) {
                retcode <- 4
              }else{
                retcode <- 7
              }
              done <- 1
            }
          }
          if( ( f<fPeak ) && (lambda>0) ) {
            fPeak <- f
            lambdaPeak <- lambda
            if( lambdaMax<=lambdaPeak ) {
              lambdaMax <- lambdaPeak*factor*factor
            }
          }
          lambda <- lambda*factor
          if( abs(lambda) > 1e20 ) {
            retcode <- 5
            done <-1
          }
        } else {
          done <- 1
          if( factor < 1.2 ){
            retcode <- 7
          } else {
            retcode <- 0
          }
        }
      }
    }
  }
  cat(sprintf("Norm of dx %10.5g\n", dxnorm))
  return(list(fhat=fhat,xhat=xhat,fcount=fcount,retcode=retcode))
}


numgrad <- function(fcn, x, ...) {
  ## fcn can return a vector, in which case numgrad returns a matrix.
  delta <- 1e-5
  ## delta <- 1e-8
  n <- length(x)
  ## we tolerate x's that may be n x 1, 1 x n, or R vectors (with no dim),
  ## but note that g comes out as n x k matrix regardless. 
  tvec <- delta*diag(n)
  f0 <- fcn(x,...)
  k <- length(f0)
  g <- matrix(0,n,k)
  badg <- FALSE
  for (i in 1:n){
    scale <- 1
    tvecv <- tvec[,i]
    if(is.null(dim(x))){
      tvecv <- as.vector(tvecv)
    }else{
      dim(tvecv) <- dim(x)
    }
    g0 <- (fcn(x+scale*tvecv,...) - f0)/(scale*delta)
    if (max(abs(g0))< 1e50){
      g[i, ] <- as.vector(g0)
    }else{
      cat("bad gradient ------------------------\n")
      badg <- TRUE
    }
  }
  return(list(g=g,badg=badg))
}


numHess <- function(fcn, x, ...) {
  f1 <- fcn
  n <- length(x)
  h <- matrix(0, n, n)
  f2 <- function(z, ...) { numgrad(fcn=f1, z, ...)$g}
  h <- numgrad(fcn=f2, x=x, ...)$g
  return(h)
}
    

SVARhtskdmdd <- function(ydata,lags,xdata=NULL, const=const, A0, lmd, Tsigbrk, breaks=NULL,
                         urprior=list(lambda=5,mu=1), mnprior=list(tight=3,decay=.5),
                         vprior=list(sig=NULL,w=1), train=0,flat=FALSE,nonorm=nonorm,ic=NULL)
### This gives the posterior integrated over A+ (the right-hand side coefficients), conditional
### on A0 and lmd.
###---------------------------------------------
### ydata:        endogenous variable data matrix, including initial condition dates.
### xdata:        exogenous variable data matrix, including initial condition dates.  
### const:        Constant term is added automatically if const=TRUE.
### A0:           Contemporaneous coefficient matrix --- constant.
### lmd:          Column Vectors of log variances of structural shocks.
### Tsigbrk:      Dates at which lmd vectors change.  Last date with old lmd (not first with new).
### breaks:       breaks in the data.  The first lags data points after a break are used
###               as new initial conditions, not data points for the fit.
### lambda:       weight on the co-persistence prior dummy observation.  (5 is reasonable)
###               lambda>0 => x variables included; lambda<0 => x variables excluded;
### mnprior       see vprior() comments
### urprior:      
### train:        If non-zero, this is the point in the sample at which the
###               "training sample" ends.  Prior x likelihood to this point is weighted to
###               integrate to 1, and therefore is treated as if it were itself the prior.
###               To do a pure training sample prior, set lambda=mu=0, mnprior=NULL, vprior$w=0,
###               train>lags.  
### flat:         Even with lambda=mu=vprior$w=0, mnprior=NULL, det(Sigma)^(-(nv+1)/2) is used
###               as a "prior", unless flat=TRUE. flat=TRUE is likely not to work unless train is reasonably large.
### nonorm:       If true, use dummy observations but do not normalize posterior to make them a
###               proper prior.  Useful to duplicate results obtained by others, to use
###               dummy observations that do not imply a proper prior, or to save computing time in case only the
###               posterior on this model's parameters, not the weight on the model, is needed.  
### ic:           Initial conditions matrix for use in forming the sums of coefficients dummy observations.
###               If ic=NULL, the means of the first lags observations in ydata are used.  If !is.null(ic),
###               ic should be a single "observation" on the y's and x's that will be used as the persistent
###               values entering the sums of coefficients dummies.
###
###               Note that to enter a prior directly as dummy observations, one can treat the
###               Dummy observations as a training sample.
###
{
    if (is.null(dim(ydata)))  ydata <- matrix(ydata, ncol=1)
    ybar <- apply(ydata[1:lags, ], 2, mean)
    T <- dim(ydata)[1]
    nv <- dim(ydata)[2]
    if (const) {
        xdata <- cbind(xdata, matrix(1,T,1))
    }
    ## looks likely that const=FALSE, xdata=NULL case crashes.  (2012.9.24)
    if (!is.null(xdata) ) stopifnot( dim(xdata)[1] == T)
    Tx <- dim(xdata)[1]
    nx <- dim(xdata)[2]
    vp <- varprior(nv,nx,lags,mnprior,vprior, urprior=urprior, ybar=ybar) # vp$: ydum,xdum,pbreaks
    ## -------- set lmd for prior dummies --------------
    if (!is.null(dim(lmd))) {
        lmdbar <- apply(lmd, 1, mean)
    } else {
        lmdbar <- lmd
    }
    ## --------------------- Tsigbrk assumed to be indexes into ydata matrix, not
    ## --------------------- dates.  Conversion from dates and adding T done in bvarWrap3().
    ## Tsigbrk <- c(invTime(Tsigbrk, ydata), T)            #dummy obs at end
    lmd <- cbind(lmd, lmdbar)
    ##-------------------------------------------
    ## var = rfvar3(ydata=rbind(ydata, vp$ydum), lags=lags, xdata=rbind(xdata,vp$xdum), breaks=matrix(c(breaks, T, T + vp$pbreaks), ncol=1),
    ## const=FALSE, lambda=lambda, mu=mu, ic=ic) # const is FALSE in this call because ones alread put into xdata
    
    

    ##Added by Joanne Im 3.5.15
     
    k=matrix(c(breaks, T, T + vp$pbreaks), ncol=1)
    
    var = rfvar3(ydata=rbind(ydata, vp$ydum), lags=lags, xdata=rbind(xdata,vp$xdum),
        breaks=k, const=FALSE, lambda=NULL,
        mu=NULL, ic=ic, sigpar=list(A0=A0,lmd=lmd,Tsigbrk=Tsigbrk))
    ##  const is FALSE in this call because ones alread put into xdata
    Tu <- dim(var$u)[1]
    if ( any(var$snglty > 0) ) error( var$snglty, " redundant columns in rhs matrix")
    lmdllh <- .5 * sum(var$lmdseries)
    llh <- -.5 * sum(var$u^2) + Tu * (-nv * log(2 * pi)/2 + determinant(A0)$modulus) +
        lmdllh
    ## nb: determinant() returns log of abs value of determinant
    nX <- lags * nv + 1
    w <-  llh + .5 * sum(var$logdetxxi) + nv * nX * log(2 * pi)/2
    if(train!=0) {
        if(train <= lags)
            {
                cat("end of training sample <= # of lags\n")  #
                    return
            }
        Tp <- train
        tbreaks <- c(breaks[breaks<train],Tp)
    } else {
        Tp <- lags
        ## because need initial conditions to form lambda/mu prior dummy obs
        tbreaks <- Tp
    }
    ytrain <- ydata[1:Tp,,drop=FALSE]
    xtrain <- xdata[1:Tp,,drop=FALSE]
    if (!nonorm) {
        priorTsigbrk <- c(0, Tp)
        ## It is assumed that there are no breaks in lmd in the training sample!
        priornsig <- 2
        priorlmd <- cbind(lmd[ , 1], lmd[ , dim(lmd)[2]])
        varp <- rfvar3(ydata=rbind(ytrain, vp$ydum), lags=lags, xdata=rbind(xtrain, vp$xdum),
                       breaks=c(tbreaks, Tp+vp$pbreaks), 
                       lambda=NULL, mu=NULL, const=FALSE, ic=ic,
                       sigpar=list(A0=A0,lmd=priorlmd, Tsigbrk=priorTsigbrk))
        ## const is FALSE here because xdata already has a column of ones.
        if (any(varp$snglty > 0)) {
            warning("Prior improper, short ", varp$snglty, " df.  Results likely nonsense.")
        } else {
            Tup <- dim(varp$u)[1]
            	#number of prior observations for the train and the dummies
            	
            lmdllhp <- .5 * sum(varp$lmdseries)
            llhp <- -.5 * sum(varp$u^2) - Tup * (nv * log(2 * pi)/2 - determinant(A0)$modulus) +
                lmdllhp
            normalizer <- .5 * sum(varp$logdetxxi) + nv * nX * log(2 * pi)/2
            wp <- llhp + normalizer
            w <- w-wp
            llh <- llh - normalizer
            ## llh is height of posterior density over A0, lmd, A+ at peak.  w is height of
            ## marginal posterior for A0, lmd, with A+ integrated out.
        }
    } else {
        varp <- NULL
    }
    return(list(w=w,var=var,varp=varp,prior=list(urprior=urprior, vprior=vprior, mnprior=mnprior)))
}


csminwelNew <- function(fcn,x0,H0,...,grad=NULL,crit=1e-7,nit,Verbose=TRUE,Long=FALSE) {
### fcn:   the objective function to be minimized.  If it has a "gradient" attribute (like output of deriv), that is used
###        as analytic gradient.  If it has a "hessian" attribute, that is used as the hessian.
### fcn0:  Lean version of fcn, without grad or hessian attributes.  May save time to provide this. (Removed this option for now.)
### x0:    initial value of the parameter vector
### H0:    initial value for the inverse Hessian.  Must be positive definite, if used.  (Not used if attr(fcn,"hessian") exists.)
### grad:  If this is a numerical vector and attr(fcn,"gradient") is not present, then grad is used as an initial gradient vector.
###        Useful for restarting if numerical gradient calculation is slow.
### crit:  Convergence criterion.  Iteration will cease when it proves impossible to improve the
###        function value by more than crit.
### nit:   Maximum number of iterations.
### ...:   A list of optional length of additional parameters that get handed off to fcn each
###        time it is called.
###        Note that if the program ends abnormally, it is possible to retrieve the current x,
###        f, and H from the files g1.mat and H.mat that are written at each iteration and at each
###        hessian update, respectively.  (When the routine hits certain kinds of difficulty, it
###        writes g2.mat and g3.mat as well.  If all were written at about the same time, any of them
###        may be a decent starting point.  One can also start from the one with best function value.)
  dots <- list(...) # (need this to save these arguments in case of cliffs)
  nx <- length(x0)
  done <- 0
  itct <- 0
  fcount <- 0
  snit <- 100
  Hlst<- list()
  badg1 <- badg2 <- badg3 <- FALSE
  f0 <- fcn(x0,...)
  NumGrad <- is.null(attr(f0,"gradient"))
  NumHess <- is.null(attr(f0,"hessian"))
  badg <- FALSE
  
  ## browser()
  if( f0 > 1e50 ) {
    stop("Bad initial parameter.")
  }
  
  ##
  #Set gradient of fcn at x0
  ##
  if( NumGrad ) {
    if( !is.numeric(grad) ) {
      gbadg <- numgrad(fcn,x0,...)
      g <- gbadg$g
      badg <- gbadg$badg
    } else {
      badg <- false
      ## This is dangerous if you use a saved g file and it
      ## turns out to have been "bad".  We used to set badg to TRUE if
      ## grad contained any zeros.
      g <- grad
    }
  } else {
    g <- attr(f0,"gradient")
    badg <- FALSE
    gbadg <- list(g=g, badg=badg)
  }
  
  retcode3 <- 101
  x <- x0
  f <- f0
  
  ##
  #Set the inverse hessian H
  ##
  if (is.null(attr(fcn,"hessian"))) {
    H <- H0
  }else{
    H <- attr(f0,"hessian")
  }
  cliff <- 0
  while( !done ) {
    g1 <- NULL; g2 <- NULL; g3 <- NULL
    ##addition fj. 7/6/94 for control
    cat('-----------------\n')
    cat('-----------------\n')
    cat(sprintf('f at the beginning of new iteration, %20.10f',f),"\n")
    if (!Long && Verbose) { # set Long=TRUE if parameter vector printouts too long
      cat("x =\n")
      print(x,digits=12)
    }
    ##-------------------------
    itct <- itct+1
    itout <- csminit(fcn,x,f,g,badg,H,...)
    f1 <- itout$fhat
    x1 <- itout$xhat
    fc <- itout$fcount
    retcode1 <- itout$retcode
    fcount <- fcount+fc
    if( retcode1 != 1 ) {         # Not gradient convergence
      ## if( retcode1==2 || retcode1==4) {
      ##   wall1 <- TRUE; badg1 <- TRUE
      ## } else {                          # not a back-forth wall, so check gradient
      if( NumGrad ) {
        gbadg <- numgrad(fcn, x1,...)
      } else {
        gbadg <- list(g=attr(f1,"gradient"),badg=FALSE)
      }
      g1 <- gbadg$g
      badg1 <- gbadg$badg
      wall1 <- (badg1 || retcode1==2 || retcode1 == 4) # A wall is back-forth line search close or bad gradient
      ## g1
      save(file="g1", g1, x1, f1, dots)
      ## }
      if( wall1 && dim(H)[1] > 1) {
        ## Bad gradient or back and forth on step length.  Possibly at
        ## cliff edge.  Try perturbing search direction, if problem is not unidimensional
        ##
        Hcliff <- H+diag(diag(H) * rnorm(nx))
        cat("Cliff.  Perturbing search direction. \n")
        itout <- csminit(fcn,x,f,g,badg,Hcliff,...)
        f2 <- itout$fhat
        x2 <- itout$xhat
        fc <- itout$fcount
        retcode2 <- itout$retcode
        fcount <- fcount+fc   # put by Jinill
        ## if(  f2 < f ) {
        ##   if( retcode2 == 2 || retcode2==4 ){
        ##     wall2 <- 1
        ##     badg2 <- 1
        ##   } else {
        if( NumGrad ){
          gbadg <- numgrad(fcn, x2,...)
        }else{
          gbadg <- list(g=attr(f2,"gradient"),badg=FALSE)
        }
        g2 <- gbadg$g
        badg2 <- gbadg$badg
        wall2 <- (badg2 || retcode2==2 || retcode2==4)
        ## g2
        ## badg2
        save(file="g2", g2, x2, f2, dots)
        if( wall2 ){
          cat('Cliff again.  Try traversing\n')
          normdx <- sqrt(sum(x2-x1)^2)
          if( normdx < 1e-13 ) { # two trial x's too close, can't traverse
            f3 <- f; x3 <- x; badg3 <- NumGrad;retcode3 <- 101
          }else{
            ## as.numeric below for robustness against f's being 1x1 arrays
            gcliff <- (x2 - x1) * (as.numeric(f2 - f1)/(normdx^2))
            dim(gcliff) <- c(nx,1)
            itout <- csminit(fcn,x,f,gcliff,0,diag(nx),...)
            f3 <- itout$fhat
            x3 <- itout$xhat
            fc <- itout$fc
            retcode3 <- itout$retcode
            fcount <- fcount+fc  # put by Jinill
            ## if( retcode3==2 || retcode3==4 ) {
            ##   wall3 <- 1
            ##   badg3 <- 1
            ## } else {
            if( NumGrad ) {
              gbadg <- numgrad(fcn, x3,...) 
            }else{
              gbadg <- list(g=attr(f3,"gradient"),badg=FALSE)
            }
            g3 <- gbadg$g
            badg3 <- gbadg$badg
            wall3 <- (badg3 || retcode3==2 || retcode3==4)
            ## g3
            ## badg3
            save(file="g3", g3, x3, f3, dots)
          }
        } else { # wall1, but not wall2, so pack f3, etc with initial values
          f3 <- f
          x3 <- x
          badg3 <- NumGrad              #i.e., use the gradient if it's analytic
          g3 <- g
          retcode3 <- 101
        }
      } else {     # no walls, or one-dimensional, so no perturbations
        f3 <- f
        x3 <- x
        badg3 <- NumGrad
        g3 <- g
        retcode3 <- 101
        f2 <- f
        x2 <- x
        badg2 <- NumGrad
        g2 <- g
        retcode2 <- 101
      }
    } else { # gradient convergence --- csminit just looked at gradient and quit
      f1 <-  f2 <-  f3 <- f; g1 <- g <- FALSE; badg2 <-  badg3 <- TRUE; retcode1 <- 1; retcode2 <- 101; retcode3 <- 101
    }
    ## how to pick gh and xh:
    ## pick first fj that improves by crit and has badg==FALSE, starting with j=3
    ## may pick other than min fj to stay away from cliff
    if( f3 < f-crit & badg3==0 ) {
      ih <- 3
      fh <- f3;xh <- x3;gh <- g3;badgh <- badg3;retcodeh <- retcode3
    } else {
      if( f2 < f-crit && badg2==0 ) {
        ih <- 2
        fh <- f2;xh <- x2;gh <- g2;badgh <- badg2;retcodeh <- retcode2
      } else {
        if( f1 < f-crit && badg1==0 ) {
          ih <- 1
          fh <- f1;xh <- x1;gh <- g1;badgh <- badg1;retcodeh <- retcode1
        } else {
          ## none qualify, so pick the one that improves most.
          fh <- min(f1,f2,f3)
          ih <- which.min(c(f1,f2,f3))
          cat("ih =",ih,"\n")
          xh <- switch(ih,x1,x2,x3)
          retcodeh <- switch(ih,retcode1,retcode2,retcode3)
          nogh <- (!exists("gh",inherits=FALSE) || is.null(gh))
          if( nogh ) {
            if( NumGrad ) {
              gbadg <- numgrad(fcn,xh,...)
            } else {
              gbadg <- list(g=switch(ih,attr(f1,"gradient"),attr(f2,"gradient"),attr(f3,"gradient")),badg=FALSE)
            }
            gh <- gbadg$g
            badgh <- gbadg$badg
          }
          badgh <- NumGrad
        }
      }
    }
    ## end of picking
    ##ih
    ##fh
    ##xh
    ##gh
    ##badgh
    stuck <- (abs(fh-f) < crit)
    if ( (!badg) && (!badgh) && (!stuck) ) {
      if(NumHess){
        H <- bfgsi(H,gh-g,xh-x)
        if(itct%%50==0){
        Hlst[[itct]] <- itout$H}
      } else {
        H <- attr(fh,"hessian")
      }
    } else {
      cat("skipped bfgsi\n")
    }
    ## if( Verbose ) {
    cat("----\n")
    cat(sprintf('Improvement on iteration %8.0f = %18.9f\n',itct,f-fh))
    ##}
    if( itct > nit ) {
      cat("iteration count termination\n")
      done <- 1
    } else {
      if( stuck ) {
        cat("improvement < crit termination\n")
        done <- 1
      }
      rc <- retcodeh
      switch( rc ,
             cat("zero gradient\n"),    #1
             cat("back and forth on step length never finished\n"), #2
             cat("smallest step still improving too slow\n"),       #3
             cat("back and forth on step length never finished\n"), #4
             cat("largest step still improving too fast\n"),        #5
             cat("smallest step still improving too slow, reversed gradient\n"), #6
             cat("warning: possible inaccuracy in H matrix\n"), #7
             )
    }
    f <- fh
    x <- xh
    g <- gh
    badg <- badgh
  }                                     # while !done
  return(list(fh=fh,xh=xh,gh=gh,H=H, Hlst=Hlst, itct=itct,fcount=fcount,retcodeh=retcodeh,...))
}
