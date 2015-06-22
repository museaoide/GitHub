
copper<-read.csv("copperfutures.csv", na.strings=c("NA"))

date01=as.Date(as.character(copper$Date01))
date03=as.Date(as.character(copper$Date03))
date06=as.Date(as.character(copper$Date06))
require(xts)
copper01=xts(copper$Settle01, date01)
copper03=xts(copper$Settle03, date03)
copper06=xts(copper$Settle06, date06)

copper01 <-na.omit(copper01)
copper03 <-na.omit(copper03)
copper06 <-na.omit(copper06)

plot.xts(copper01)
plot.xts(copper03)
plot.xts(copper06)

copper0301=copper03-copper01
copper0603=copper06-copper03

pdf("copper0301.pdf")
plot.xts(copper0301, type="h")
abline(h=0)
dev.off()


pdf("copper0603.pdf")
plot.xts(copper0603, type="h")
abline(h=0)
dev.off()

gdp<-read.csv("gdp.csv", na.strings=c("NA"))