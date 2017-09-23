# Final copy
# pwr_comparison_1
library(gtools)


data1= data.frame()
Nobel<-function(n){
  p11hat<-rep(0,n)
  p12hat<-rep(0,n)
  p21hat<-rep(0,n)
  p22hat<-rep(0,n)
  Rho<- rep(0,n)
  AsyVarRho<-rep(0,n)
  OddsRatio<-rep(0,n)
  lnOR<-rep(0,n)
  AsyVarOR<-rep(0,n)
  count_lnOR<-rep(0,n)
  count_Rho<-rep(0,n)


  for (i in 1:n) {
    MB <- rdirichlet(1, c(1,1,1,1))
    Sample<-rmultinom(n, 100, prob= c(MB[1,1], MB[1,2], MB[1,3], MB[1,4]))

    p11hat <- Sample[1,i]/100
    p12hat <- Sample[2,i]/100
    p21hat <- Sample[3,i]/100
    p22hat <- Sample[4,i]/100


    p1p <- p11hat + p12hat
    p2p <- 1-p11hat-p12hat
    pp1 <- p11hat+p21hat
    pp2 <- 1-p11hat-p21hat
    Var1 <- p11hat*(1-p11hat)
    Var2 <- p12hat*(1-p12hat)
    Var3 <- p21hat*(1-p21hat)
    Cov12 <- -p11hat*p12hat
    Cov13 <- -p11hat*p21hat
    Cov23 <- -p12hat*p21hat
    U <- sqrt(p1p*p2p*pp1*pp2)
    V <- p11hat*p2p*pp2 - p12hat*p2p*pp1 - p21hat*p1p*pp2 + (1-p11hat-p12hat-p21hat)*p1p*pp1

    DUDx <- (1-p1p-pp1)*(p1p+pp1-2*p1p*pp1)
    DUDy <- pp1*(1-pp1)*(1-2*p1p)
    DUDz <- p1p*(1-p1p)*(1-2*pp1)
    DVDx <- 1-p1p-pp1
    DVDy <- -pp1
    DVDz <- -p1p
    DfDx <- (1/sqrt(U))*DVDx - 0.5*U^(-1.5)*V*DUDx
    DfDy <- (1/sqrt(U))*DVDy - 0.5*U^(-1.5)*V*DUDy
    DfDz <- (1/sqrt(U))*DVDz - 0.5*U^(-1.5)*V*DUDz
    AsyVarRho <- Var1*DfDx^2 + Var2*DfDy^2 + Var3*DfDz^2 + 2*Cov12*DfDx*DfDy +2*Cov13*DfDx*DfDz + 2*Cov23*DfDy*DfDz

    OddsRatio <- p11hat*(1-p11hat-p12hat-p21hat)/(p12hat*p21hat)
    lnOR <- log(OddsRatio)

    AsyVarOR <- (1/p11hat + 1/p12hat + 1/p21hat + 1/(1 - p11hat - p12hat - p21hat))
    Rho <- V/sqrt(U)
    z1<-lnOR/sqrt(AsyVarOR/100)
    z2<- Rho/ sqrt(AsyVarRho/100)
    count_lnOR<-ifelse(abs(z1)>1.96,1,0)
    count_Rho<-ifelse(abs(z2)>1.96,1,0)

    df <- c(p11hat, p12hat, p21hat, Rho, AsyVarRho, OddsRatio,lnOR, AsyVarOR,count_lnOR,count_Rho)
    data1=rbind(df,data1)

  }# for loop

  colnames(data1)<-c("p11hat", "p12hat", "p21hat", "Rho", "AsyVarRho", "OddsRatio","lnOR","AsyVarOR","count_lnOR",
                     "count_Rho")


  return(data1)
}

Y<-Nobel(10)
Y



C_lnOR<-as.matrix(table( Y$count_lnOR))
C_lnOR

C_Rho<-as.matrix(table( Y$count_Rho))
C_Rho
Percent_lnOR<- (C_lnOR[2]/ sum(C_lnOR[1], C_lnOR[2]))*100
Percent_lnOR
Percent_Rho<- (C_Rho[2]/ sum(C_Rho[1], C_Rho[2]))*100
Percent_Rho


result<-replicate(n = 10, expr = {
  Y <- Nobel(1000)
  C_lnOR <- as.matrix(table(Y$count_lnOR))
  C_Rho <- as.matrix(table(Y$count_Rho))
  Percent_lnOR <- (C_lnOR[2] / sum(C_lnOR[1], C_lnOR[2])) * 100
  Percent_Rho <- (C_Rho[2] / sum(C_Rho[1], C_Rho[2])) * 100
  c(Percent_lnOR, Percent_Rho)
})

row.names(result)<- c("Power_lnOR", "Power_Rho")
result
