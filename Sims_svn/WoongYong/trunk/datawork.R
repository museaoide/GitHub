load("dvardata2.RData")
load("USDa.Rdata")
load("FedRcptExpT.RData")
load("FFRates.RData")

PriSurp1 <- FedRcptExpT[,46] - (FedRcptExpT[,13] - FedRcptExpT[,28])
PriSurp2 <- FedRcptExpT[,37] - rowSums(FedRcptExpT[,41:44],na.rm=TRUE) + FedRcptExpT[,45] + FedRcptExpT[,28]

PriSurp <- ts(c(window(PriSurp2,start=c(1947,1),end=c(1959,4)),window(PriSurp1,start=c(1960,1))),start=c(1947,1),frequency=4)
tau <- PriSurp / exp(dvardata2[,"Cons.Defl"]) * 100

y <- window(cbind(r=FF/4,p=dvardata2[,"Cons.Defl"],c=dvardata2[,"GDP.Real"],tau=tau),start=c(1954,4),end=c(2005,1))

b <- USDa[,"MktblTrsyDbtMkt"] / exp(dvardata2[,"Cons.Defl"]) * 100

y2 <- window(cbind(r=FF/4,p=dvardata2[,"Cons.Defl"],c=dvardata2[,"GDP.Real"],tau=tau,b=b),start=c(1954,4),end=c(2005,1))

# y0 <- window(cbind(r=FF/4,p=dvardata2[,"Cons.Defl"],c=dvardata2[,"GDP.Real"],tau=tau),start=c(1954,3),end=c(1954,3))
# b0 <- window(USDa[,"MktblTrsyDbtMkt"],start=c(1954,3),end=c(1954,3)) / exp(y0[,"p"]) * 100

#                 r        p        c       tau
#1954 Q3 0.002566667 2.920847 7.633788 -1.616641
# b0 830.9533

# y0 <- window(cbind(r=FF/4,p=dvardata2[,"Cons.Defl"],c=dvardata2[,"GDP.Real"],tau=tau),start=c(1982,4),end=c(1982,4))
# b0 <- window(USDa[,"MktblTrsyDbtMkt"],start=c(1982,4),end=c(1982,4)) / exp(y0[,"p"]) * 100

#                 r        p        c       tau
#1982 Q4 0.02321667 4.110415 8.554458 -150.2329
# b0 1412.998

# FFMonth <- ts(read.table(file="FFRates.dat"),start=c(1954,7),frequency=12)
# FF <- ts(t(kronecker(diag(nrow=length(FFMonth)/3),c(1,1,1))) %*% FFMonth / 3,start=c(1954,3),frequency=4)