screpD <- function(pgw=rep(0,times=nx*ny-nx*(nx-1)/2-ny),nx=8,ny=16,U=function(a,b){log(ifelse(a <= 0,1e-100,a))+log(ifelse(b-a <= 0,1e-100,b-a))},xrange=c(0,1),wrange=c(0,1),gf=function(z)sin(z/pi),lambda=1)
  {
    ## computes objective function E[U] minus lambda * (information constraint)
    ## pgw:   defines pdf of x given w. The pdf is defined on a matrix of x and w values, but passed as a vector
    ##        to this function.  pdf matrix eventually has nx rows, but only nx-1 are in pgw, due to
    ##        normalization. Since we consider only U's defined on x<=w, lower triangle of pgw is omitted. 
    ## nx:    number of discrete x values considered (and rows in pgw)
    ## U:     objective function U(x,w)
    ## xrange: range of values for x (grid is found from this and nx)
    ## wrange: range of values for w (grid is found from this and length(pgw)/nx
    ## gg:     marginal pdf of w
    ## lambda: multiplier on the information constraint, cost of capacity
    ##
    lbp <- 1e-7
    pnlty <- 1
    x <- seq(xrange[1],xrange[2],length.out=nx)
    w <- seq(wrange[1],wrange[2],length.out=ny)
    gg <- gf(x)
    cltw <- outer(x,w,function(a,b){a<b})
    pgwLoc <- matrix(-100,nx-1,ny)
    pgwLoc[cltw[2:nx,]] <- pgw
    ##-------------------------------------------------------------------------------
    ## Map pgw into a matrix whose values are all in (0,1) and sum to one in each col.
    pgwLoc <- exp(pgwLoc)
    pgwLoc <- rbind(rep(1,ny), pgwLoc)
    colsumi <- 1/(as.vector(rep(1,nx) %*% pgwLoc))
    pgwLoc <- t(colsumi * t(pgwLoc))
    ##--------------------------------------------------------------------------------
    pdf <- t(gg * t(pgwLoc))
    EU <- sum(pdf * outer(x,w,FUN=U))
    pc <- pdf %*% rep(1,dim(pdf)[2])
    k <- sum(pdf * log(ifelse(pdf>0,pdf,1))) - sum(pc * log(ifelse(pc>0,pc,1))) - sum(gg * log(ifelse(gg>0,gg,1)))
    ## obf <- EU - lambda * k
    ## two lines below apply a penalty for pdf values near zero.  This prevents the search from wandering over large
    ## changes in large negative values of log p's, and may speed convergence.
    ##---------------------------------
    browser()
    zerop <- pdf[cltw & (pdf < lbp)]
    obf <- EU - lambda *k - sum(pnlty * (lbp / zerop - 1))
    ##---------------------------------
    return(list(obf=obf,EU=EU,pc=pc,pdf=pdf,pgw=pgwLoc,info=k/log(2)))
  }
