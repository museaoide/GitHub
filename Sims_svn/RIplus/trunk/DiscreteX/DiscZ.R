DiscZ <- function(x, nx, gy, y, U, xfcn, alph, ...) {
    ##
    ## xfcn:   A function of ygivenx (cond'l prob matrix), y and U returning optimal x
    ## There is a special version for CapM, because need to impose that rows of x sum to one in that model
    ## Must be able to evaluate with x a matrix, with each col an arg
    require(abind)
    if (is.null(dim(x)) ) {
        nargcol <- 1
        x <- matrix(x, ncol=1)
    } else {
        nargcol <- dim(x)[2]
    }
    pold <- matrix(x[1:(nx-1), ], ncol=nargcol)
    pold <- rbind(pold, 1 - apply(pold, 2, sum))
    ny <- dim(y)[2]
    ## if capm
    ## xold <- array(x[-(1:(nx-1)), ], c(nx, ny - 1, nargcol ))
    ## sxold <-  1 - apply(xold, c(1,3), sum)
    ## sxold <- array(sxold, c(nx , 1, nargcol))
    ## xold <- abind(xold, sxold, along=2)
    ## else
    xold <- array(x[nx:(nx * ny - 1), ], c(nx, ny, nargcol))
    pnew <- matrix(0, nx, nargcol)
    xnew <- array(0, c(nx, ny, nargcol))
    screp <- matrix(1e20, dim(x)[1], nargcol)
    ## DiscPObjXmv needs pold truncated, but not xold
    for (icol in 1:nargcol) {
        if (!any(pold[ , icol] < 0)) {
            xc <- c(nx, pold[1:(nx - 1), icol], xold[ , , icol])
            Dout <- DiscPObjXmv(xc, gy, y, U=U, alph=alph,...)    
            pnew[ , icol] <- Dout$pnew
            ## xnew <- Dout$ygivenx %*% y #-------xnew solves E[Dxu %*% y = 0] This is for UmvTrack
            xnew[ , , icol] <- xfcn(Dout$ygivenx, y)
            ## if !capm
            screp[ , icol] <- c((pnew[-nx , icol] - pold[-nx , icol]), xnew[ , , icol] - xold[ , , icol])
            ## else
            ## screp[ , icol] <- c((pnew[-nx , icol] - pold[-nx , icol]), (xnew - xold)[ ,1:(ny-1), icol] )
            ## fi
            ##otherwise leave the 1e20 values in place in that column
        }
    }
    ## x <- c(nx, pnew[-nx], xnew)
    ## if (itct %% 100 == 0) {
    ##   print(itct)
    ##   print(screp)
    ##   print(Dout$obj)
    ##   print(Dout$pnew)
    ##   print(xnew)
    ## }
    return(screp)
}

