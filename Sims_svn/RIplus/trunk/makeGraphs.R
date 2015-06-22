makeGraphs <- function(csout,U,lambda,gamma,gg, tag="") {
  scrout <- screpD(csout$xh, nx=16, ny=31, U=U, xrange=c(.01,.5), wrange=c(.02,1), gg=gg, lambda=lambda)
  with(
       scrout, {
         gtk()
         plot(x2, pc, type="l", main=bquote(paste("pdf of consumption, ", .(tag), ~~~~gamma==.(gamma)~~~~ lambda==.(lambda),sep=" ")), xlab="C", ylab="p(c)")
         dev.copy2eps(file=paste("g",as.character(gamma),"l",as.character(lambda),tag,"pc.eps",sep=""))
         pgc <- as.vector(1/pc) * pdf
         gtk(); image(x2,w2,pgw,col=gray(64:1/64),main=bquote(paste("p(c|w), ", .(tag), ~~~~gamma==.(gamma)~~~~ lambda==.(lambda))), xlab="C",ylab="W")
         lines(x2,x2,lty="dashed"); lines(x2,x2*2)
         dev.copy2eps(file=paste("g",as.character(gamma),"l",as.character(lambda),tag,"pgw.eps",sep=""))
         gtk();  image(x2,w2,pdf,col=gray(64:1/64),main=bquote(paste("p(c,w), ", .(tag), ~~~~gamma==.(gamma)~~~~ lambda==.(lambda))), xlab="C",ylab="W")
         lines(x2,x2,lty="dashed"); lines(x2,x2*2)
         dev.copy2eps(file=paste("g",as.character(gamma),"l",as.character(lambda),tag,"pdf.eps",sep=""))
         gtk();  image(x2,w2,pgc,col=gray(64:1/64),main=bquote(paste("p(w|c), ", .(tag), ~~~~gamma==.(gamma)~~~~ lambda==.(lambda))), xlab="C",ylab="W")
         lines(x2,x2,lty="dashed"); lines(x2,x2*2)
         dev.copy2eps(file=paste("g",as.character(gamma),"l",as.character(lambda),tag,"pgc.eps",sep=""))
       }
       )
  return(scrout)
}
