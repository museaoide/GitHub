plotir <- function(resp, var, shock, type=c("multiple","single"),main="Impulse response plot", yax.flip=TRUE, vfirst=TRUE) {
    ## resp: impulse responses, as produced by impulsdtrf
    ## var: indexes of variables whose responses are to be plotted
    ## shock:  indexes of shocks whose effects are to be plotted
    ## yax.flip: Unless very few responses are plotted, vertical axis
    ##           labels start to overlap, unless yax.flip=TRUE
    ## vfirst: the variable index increases down each column.  Otherwise
    ##         (vfirst=FALSE) the shock index does so.
    ## As opposed to more direct plotting commands, this one has the
    ## advantage that (1) it scales all responses of the same variable
    ## to the same height on a given graph, (2) that it allows printing
    ## a matrix of graphs with rows corresponding to variables and columns
    ## to shocks, or vice versa, and (3) that it prompts for the legend
    ## (in the "single" case).
    ##-------------------------------
    ## Program seems to be unfinished or truncated somehow.  Incomplete repair 13.10.26.
    pr <- resp[var,shock, , drop=FALSE]
    pdf(file="irplot.pdf", width=21, height=24)
    ybound <- matrix(0,2,length(var))
    for (iv in 1:length(var)) ybound[ , iv] <- range(pr[iv, , ])
    if(!vfirst) pr <- aperm(pr, c(2,1,3))
    nr <- dim(pr)[1]
    nc <- dim(pr)[2]
    np <- dim(pr)[3]
    par(omi=c(1,1,2,1))
    op <- par(mfrow=c(nr,nc))
    for (ir in 1:nr) {
        for (ic in 1:nc) {
            if (ic == 1) {
                par(mai=c(1/nr,6/nc,2/nr,1/nc))
            } else {
                par(mai=c(1/nr,2/nr,2/nr,1/nc))
            }
            plot(0:(np-1), pr[ir, ic, ],bty="l", type="l",
                 ylim=ybound[ , if(vfirst) ir else ic], frame=TRUE,
                 ylab=if(ic==1) dimnames(pr)[[1]][ir] else "", xlab="",
                 main=if(ir ==1) dimnames(pr)[[2]][ic] else "",
                 xaxt=if(ir==nr)"s" else "n",
                 yaxt=if(ic==1)"s" else "n", cex.main=2, cex.lab=2)
            lines(c(0,(np-1)),c(0,0)) #x axis line
        }
    }
    title(main, outer=TRUE, cex.main=2)        
    dev.off()
    par(op)
}

