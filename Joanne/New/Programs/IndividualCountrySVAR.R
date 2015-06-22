
listData6e3_aus=tvvData_noffr4(data=data, commdata=commdata,commIndexWeights=c(0,0,1,0), countries=c("aus"), begin=NULL, terminate=NULL, SecondaryBrk=c(0, 1984, 1998.75, 2006.75, 2015))
results6e2_aus=plutusTvv3(listData6e3_aus,nit=600,H1=10^-2,H2=10^-3)
save(results6e2_aus, file="results6e2_aus.RData")

listData6e3_be=tvvData_noffr4(data=data, commdata=commdata,commIndexWeights=c(0,0,1,0), countries=c("be", "be_posteuro"), begin=NULL, terminate=NULL, SecondaryBrk=c(0, 1984, 1998.75, 2006.75, 2015))
results6e2_be=plutusTvv3(listData6e3_be,nit=600,H1=10^-2,H2=10^-3)
save(results6e2_be, file="results6e2_be.RData")

listData6e3_can=tvvData_noffr4(data=data, commdata=commdata,commIndexWeights=c(0,0,1,0), countries=c("can"), begin=NULL, terminate=NULL, SecondaryBrk=c(0, 1984, 1998.75, 2006.75, 2015))
results6e2_can=plutusTvv3(listData6e3_can,nit=600,H1=10^-2,H2=10^-3)
save(results6e2_can, file="results6e2_can.RData")

listData6e3_che=tvvData_noffr4(data=data, commdata=commdata,commIndexWeights=c(0,0,1,0), countries=c("che"), begin=NULL, terminate=NULL, SecondaryBrk=c(0, 1984, 1998.75, 2006.75, 2015))
results6e2_che=plutusTvv3(listData6e3_che,nit=600,H1=10^-2,H2=10^-3)
save(results6e2_che, file="results6e2_che.RData")

listData6e3_dan=tvvData_noffr4(data=data, commdata=commdata,commIndexWeights=c(0,0,1,0), countries=c("dan"), begin=NULL, terminate=NULL, SecondaryBrk=c(0, 1984, 1998.75, 2006.75, 2015))
results6e2_dan=plutusTvv3(listData6e3_dan,nit=600,H1=10^-2,H2=10^-3)
save(results6e2_dan, file="results6e2_dan.RData")

listData6e3_fr=tvvData_noffr4(data=data, commdata=commdata,commIndexWeights=c(0,0,1,0), countries=c("fr", "fr_posteuro"), begin=NULL, terminate=NULL, SecondaryBrk=c(0, 1984, 1998.75, 2006.75, 2015))
results6e2_fr=plutusTvv3(listData6e3_fr,nit=600,H1=10^-2,H2=10^-3)
save(results6e2_fr, file="results6e2_fr.RData")

listData6e3_it=tvvData_noffr4(data=data, commdata=commdata,commIndexWeights=c(0,0,1,0), countries=c("it", "it_posteuro"), begin=NULL, terminate=NULL, SecondaryBrk=c(0, 1984, 1998.75,2006.75,2015))
results6e2_it=plutusTvv3(listData6e3_it,nit=600,H1=10^-3,H2=10^-4)
save(results6e2_it, file="results6e2_it.RData")

listData6e3_jpn=tvvData_noffr4(data=data, commdata=commdata,commIndexWeights=c(0,0,1,0), countries=c("jpn"), begin=NULL, terminate=NULL, SecondaryBrk=c(0, 1984, 1998.75,2006.75,2015))
results6e2_jpn=plutusTvv3(listData6e3_jpn,nit=600,H1=10^-2,H2=10^-3)
save(results6e2_jpn, file="results6e2_jpn.RData")

listData6e3_nl=tvvData_noffr4(data=data, commdata=commdata,commIndexWeights=c(0,0,1,0), countries=c("nl", "nl_posteuro"), begin=NULL, terminate=NULL, SecondaryBrk=c(0, 1984, 1998.75,2006.75,2015))
results6e2_nl=plutusTvv3(listData6e3_nl,nit=600,H1=10^-2,H2=10^-3)
save(results6e2_nl, file="results6e2_nl.RData")

listData6e3_no=tvvData_noffr4(data=data, commdata=commdata,commIndexWeights=c(0,0,1,0), countries=c("no"), begin=NULL, terminate=NULL, SecondaryBrk=c(0, 1984, 1998.75, 2006.75, 2015))
results6e2_no=plutusTvv3(listData6e3_no,nit=600,H1=10^-2,H2=10^-3)
save(results6e2_no, file="results6e2_no.RData")

listData6e3_pt=tvvData_noffr4(data=data, commdata=commdata,commIndexWeights=c(0,0,1,0), countries=c("pt", "pt_posteuro"), begin=NULL, terminate=NULL, SecondaryBrk=c(0, 1984, 1998.75, 2006.75, 2015))
results6e2_pt=plutusTvv3(listData6e3_pt,nit=600,H1=10^-3,H2=10^-5)
save(results6e2_pt, file="results6e2_pt.RData")

listData6e3_se=tvvData_noffr4(data=data, commdata=commdata,commIndexWeights=c(0,0,1,0), countries=c("se"), begin=NULL, terminate=NULL, SecondaryBrk=c(0, 1984, 1998.75, 2006.75, 2015))
results6e2_se=plutusTvv3(listData6e3_se,nit=600,H1=10^-2,H2=10^-3)
save(results6e2_se, file="results6e2_se.RData")

listData6e3_uk=tvvData_noffr4(data=data, commdata=commdata,commIndexWeights=c(0,0,1,0), countries=c("uk"), begin=NULL, terminate=NULL, SecondaryBrk=c(0, 1984, 1998.75, 2006.75, 2015))
results6e2_uk=plutusTvv3(listData6e3_uk,nit=600,H1=10^-2,H2=10^-3)
save(results6e2_uk, file="results6e2_uk.RData")


listData6e3_us=tvvData_noffr4(data=data, commdata=commdata,commIndexWeights=c(0,0,1,0), countries=c("us"), begin=NULL, terminate=NULL, SecondaryBrk=c(0, 1984, 1998.75,2006.75,2015))
results6e2_us=plutusTvv3(listData6e3_nl,nit=600,H1=10^-2,H2=10^-3)
save(results6e2_us, file="results6e2_us.RData")


listData6e3_za=tvvData_noffr4(data=data, commdata=commdata,commIndexWeights=c(0,0,1,0), countries=c("za"), begin=NULL, terminate=NULL, SecondaryBrk=c(0, 1984, 1998.75, 2006.75, 2015))
results6e2_za=plutusTvv3(listData6e3_za,nit=600,H1=10^-2,H2=10^-3)
save(results6e2_za, file="results6e2_za.RData")



IPaus=IPTvv(results6e2_aus, errorbands=FALSE, nit=20000)
IPbe=IPTvv(results6e2_be, errorbands=FALSE, nit=20000)
IPcan=IPTvv(results6e2_can, errorbands=FALSE, nit=20000)
IPche=IPTvv(results6e2_che, errorbands=FALSE, nit=20000)
IPdan=IPTvv(results6e2_dan, errorbands=FALSE, nit=20000)
IPfr=IPTvv(results6e2_fr, errorbands=FALSE, nit=20000)
IPit=IPTvv(results6e2_it, errorbands=FALSE, nit=20000)
IPjpn=IPTvv(results6e2_jpn, errorbands=FALSE, nit=20000)
IPnl=IPTvv(results6e2_nl, errorbands=FALSE, nit=20000)
IPno=IPTvv(results6e2_no, errorbands=FALSE, nit=20000)
IPpt=IPTvv(results6e2_pt, errorbands=FALSE, nit=20000)
IPse=IPTvv(results6e2_se, errorbands=FALSE, nit=20000)
IPuk=IPTvv(results6e2_uk, errorbands=FALSE, nit=20000)
IPus=IPTvv(results6e2_us, errorbands=FALSE, nit=20000)
IPza=IPTvv(results6e2_za, errorbands=FALSE, nit=20000)

