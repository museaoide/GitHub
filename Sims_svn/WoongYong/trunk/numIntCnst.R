numIntCnst <- function(A,m=30,delta=1){
# numerically calculates the integral of a matrix exponential. 
# \int_0^delta exp(-As) ds
#
# delta: the end point of the integration
# m    : the interval [0,delta] will be divided into 2^(m-1) subintervals
#        (keep 1<m<40)

mint  <- delta*2^(1-m);
expminA <- padm(-A*mint)
intgr <- expminA
for (indx in 2:m){
	intgr <- intgr + intgr%*%expminA
	expminA <- expminA%*%expminA
}
intgr <- mint*intgr
return(intgr)
} #---------------------------------------- THE END