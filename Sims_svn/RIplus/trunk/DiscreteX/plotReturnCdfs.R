plotReturnCdfs <- function(Dout, alph, ptfno, y) {
    xlim <- range(y)
    ybar <- mean(y[ , -1]) # first col is constant, others assumed drawn from same distribution
    ysd <- sd(y[ , -1])
    pn <- pnorm(seq(xlim[1], xlim[2], length=1000), mean=ybar, sd=ysd)
    calph <- as.character(alph)
    cptfno <- as.character(ptfno)
    plot(seq(xlim[1], xlim[2], length=1000),pn, main=substitute(paste("cdf's for total returns, ",alpha, " = ", alph, ", portfolio ", ptfno), list(alph=alph, ptfno=ptfno)), xlab="total return", ylab="probability", typ="l")
    lines(cdfwtd(Dout$ygivenx[ptfno, ], y[ , 2]), col=2)
    lines(cdfwtd(Dout$ygivenx[ptfno, ], y[ , 3]), col=3)
    grid()
    legend("topleft", legend=c("unconditional", "asset 2", "asset 3"), col=1:3, lwd=2)
}
         
