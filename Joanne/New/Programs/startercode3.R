

rfvar4 <- function(listData=NA,  lcIndA0=lcIndA0, lcIntA0=lcIntA0, ydata=NA,lags=nLags,xdata=NULL,const=const,breaks=NULL,lambda=5,mu=2,ic=NULL, sigpar=NULL, x) {
    ## This algorithm goes for accuracy without worrying about memory requirements.
    ## ---------------------------------------------------------------------------
    ## The standard prior it implements is NOT APPROPRIATE for seasonally unadjusted data, even
    ## if seasonal dummies are included in xdata.  The prior shrinks toward simple persistence, so it
    ## will tend to prevent the dummies from picking up all the seasonality.
    ## ---------------------------------------------------------------------------
    ##listData: output from tvvData_* programs 
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
    ## sigpar: list(lmd, Tsigbrk) Allow SVAR with time varying shock variances. Tsigbrk, if not null, must begin with 0. See below.
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
        
        
        ## Allow for country specific A0
       
		nCountries=length(unique(listData$countries))
		    
		    
		CountryRegimes=cumsum(listData$dlength0-lags)
		T=dim(y)[1]
		CountryRegimes =c(0, CountryRegimes, T)
		
		
	
		nIntA<-sum(lcIntA0)
		nIndA<-sum(lcIndA0)
		nVar<-dim(listData$Y)[2]
		

        
        ya0=NULL
        
        #+1 because we have our priors which are appended to the end
        
        for (i in 1:(nCountries+1)){
    		iIntA0=matrix(0, nVar, nVar)
        	iIndA0=matrix(0, nVar, nVar)
        	iIntA0[lcIntA0]=x[1:nIntA]
   		iIndA0[lcIndA0]=x[nIntA+(i-1)*nIndA + (1: nIndA)]
   		iA0=iIntA0+iIndA0
   		diag(iA0)=1
   		
   		if (i==nCountries+1){
   			iA0=iIntA0
   			diag(iA0)=1
   		}	
   			
        	iy=y[(CountryRegimes[i]+1): CountryRegimes[i+1],]
        	ya0 <- rbind(ya0, iy %*% t(iA0))
        	
        }
        
     

        B <- matrix(0,  nX, nvar)
        u <- matrix(0, Tsmpl, nvar)
        uraw <- u
        xxi <- array(0, c(nX, nX, nvar))
        logdetxxi <- vector("numeric", nvar)
        snglty <- vector("numeric", nvar)
        for (iq in 1:nvar) {
            wt <- exp(.5 * lmdseries[iq, ])
   #note that the * isn't matrix multiplication
            Xq <-  wt * X
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



SVARhtskdmdd3 <- function(listData, lcIndA0=lcIndA0, lcIntA0=lcIntA0, ydata,lags,xdata=NULL, const=const, iAbar, x, lmd=lmd, Tsigbrk, Countrybrk, breaks=NULL,
                         urprior=list(lambda=5,mu=1), mnprior=list(tight=3,decay=.5),
                         vprior=list(sig=NULL,w=1),  train=0,flat=FALSE,nonorm=nonorm,ic=NULL)
### Function: This gives the posterior integrated over A+ (the right-hand side coefficients), conditional
### on A0 and lmd. Modified for plutusTvv3, bvarwrap6....
### Modified such that priors are implemented for each country

###---------------------------------------------
### ydata:        endogenous variable data matrix, including initial condition dates.
### xdata:        exogenous variable data matrix, including initial condition dates.  
### const:        Constant term is added automatically if const=TRUE.
### A0:           Contemporaneous coefficient matrix --- constant.
### lmd:          Column Vectors of log variances of structural shocks.
### Tsigbrk:      Dates at which lmd vectors change.  Last date with old lmd.
### Countrybrk:   Dates before which series is of another country. 
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
    xbar <- apply(xdata, 2, mean)
    #MARGIN: a vector giving the subscripts which the function will be applied over. E.g., for a matrix 1 indicates rows, 2 indicates columns, c(1, 2) indicates rows and columns. Where X has named dimnames, it can be a character vector selecting dimension names.
    
    T <- dim(ydata)[1]
    nv <- dim(ydata)[2]
    
    if (const) {
        xdata <- cbind(xdata, matrix(1,T,1))
    }
    ## looks likely that const=FALSE, xdata=NULL case crashes.  (2012.9.24)
    if (!is.null(xdata) ) stopifnot( dim(xdata)[1] == T)
    Tx <- dim(xdata)[1]
    nx <- dim(xdata)[2]
    vp <- varprior(nv,nx,lags,mnprior,vprior, urprior=urprior, ybar=ybar, xbar=xbar) # vp$: ydum,xdum,pbreaks
    ## -------- form lmd --------------
  	
	
   
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
    breaks=Countrybrk
    
    k=matrix(c(breaks, T, T + vp$pbreaks), ncol=1)
    var = rfvar4(listData=listData, lcIndA0= lcIndA0, lcIntA0= lcIntA0, ydata=rbind(ydata, vp$ydum), lags=lags, xdata=rbind(xdata,vp$xdum),
        breaks=k, const=FALSE, lambda=NULL, mu=NULL, ic=ic, x=x, sigpar=list(lmd=lmd,Tsigbrk=Tsigbrk))
        
      
  
    ##  const is FALSE in this call because ones alread put into xdata
    Tu <- dim(var$u)[1]
    if ( any(var$snglty > 0) ) error( var$snglty, " redundant columns in rhs matrix")
    lmdllh <- .5 * sum(var$lmdseries)
    llh <- -.5 * sum(var$u^2) + Tu * (-nv * log(2 * pi)/2 + determinant(iAbar)$modulus) +
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
    
    ##We have nonorm=
    if (!nonorm) {
    	 
        priorTsigbrk <- c(0, Tp)
        ## It is assumed that there are no breaks in lmd in the training sample!
        priornsig <- 2
        priorlmd <- cbind(lmd[ , 1], lmd[ , dim(lmd)[2]])
        varp <- rfvar3(ydata=rbind(ytrain, vp$ydum), lags=lags, xdata=rbind(xtrain, vp$xdum),
                       breaks=c(tbreaks, Tp+vp$pbreaks), 
                       lambda=NULL, mu=NULL, const=FALSE, ic=ic,
                       sigpar=list(A0=iAbar,lmd=priorlmd, Tsigbrk=priorTsigbrk))
        ## const is FALSE here because xdata already has a column of ones.
        if (any(varp$snglty > 0)) {
            warning("Prior improper, short ", varp$snglty, " df.  Results likely nonsense.")
        } else {
        	
        	
            Tup <- dim(varp$u)[1]
            	#number of prior observations for the train and the dummies
            	
            #lmd prior
            lmdllhp <- .5 * sum(varp$lmdseries)
            	#remember that lmd is -log(variance)
            
            #prior
         
            llhp <- -.5 * sum(varp$u^2) - Tup * (nv * log(2 * pi)/2 - determinant(iAbar)$modulus) +
                lmdllhp 
                                
            normalizer <- .5 * sum(varp$logdetxxi) + nv * nX * log(2 * pi)/2
            
            wp <- llhp + normalizer
           
           
           #A0, lmd, with A+ integrated out.
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

