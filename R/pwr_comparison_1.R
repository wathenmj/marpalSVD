# Final copy
# pwr_comparison_1
library(gtools)
# create output dateframe
cell_counts_data= data.frame()
power_function<-function(n){
  # generate joint distribution estimates
  p11hat<-rep(0,n)
  p12hat<-rep(0,n)
  p21hat<-rep(0,n)
  p22hat<-rep(0,n)
  Rho<- rep(0,n);
  AsyVarRho<-rep(0,n) #canonical
  OddsRatio<-rep(0,n)
  lnOR<-rep(0,n)
  AsyVarOR<-rep(0,n) #odds ratio
  count_lnOR<-rep(0,n)
  count_Rho<-rep(0,n)


  for (i in 1:n) {
    # generating the joint distribution
    MB <- rdirichlet(1, c(1,1,1,1))
    cell_counts<-rmultinom(n, 100, prob= c(MB[1,1], MB[1,2], MB[1,3], MB[1,4]))
    # generate joint distribution estimates
    p11hat <- cell_counts[1,i]/100
    p12hat <- cell_counts[2,i]/100
    p21hat <- cell_counts[3,i]/100
    p22hat <- cell_counts[4,i]/100
    # row and column marginals
    p1p <- p11hat + p12hat
    p2p <- 1-p11hat-p12hat
    pp1 <- p11hat+p21hat
    pp2 <- 1-p11hat-p21hat
    # variance and covariance of estimates
    Var1 <- p11hat*(1-p11hat)
    Var2 <- p12hat*(1-p12hat)
    Var3 <- p21hat*(1-p21hat)
    Cov12 <- -p11hat*p12hat
    Cov13 <- -p11hat*p21hat
    Cov23 <- -p12hat*p21hat
    U <- sqrt(p1p*p2p*pp1*pp2) # for simplification of denominator and numerator
    V <- p11hat*p2p*pp2 - p12hat*p2p*pp1 - p21hat*p1p*pp2 + (1-p11hat-p12hat-p21hat)*p1p*pp1
    # generating canonical correlation Rho
    Rho <- V/sqrt(U)
    # partial derivatives
    DUDx <- (1-p1p-pp1)*(p1p+pp1-2*p1p*pp1)
    DUDy <- pp1*(1-pp1)*(1-2*p1p)
    DUDz <- p1p*(1-p1p)*(1-2*pp1)
    DVDx <- 1-p1p-pp1
    DVDy <- -pp1
    DVDz <- -p1p
    DfDx <- (1/sqrt(U))*DVDx - 0.5*U^(-1.5)*V*DUDx
    DfDy <- (1/sqrt(U))*DVDy - 0.5*U^(-1.5)*V*DUDy
    DfDz <- (1/sqrt(U))*DVDz - 0.5*U^(-1.5)*V*DUDz
    # asymptotic variance of canonical correlation Rho
    AsyVarRho <- Var1*DfDx^2 + Var2*DfDy^2 + Var3*DfDz^2 + 2*Cov12*DfDx*DfDy +2*Cov13*DfDx*DfDz + 2*Cov23*DfDy*DfDz

    OddsRatio <- p11hat*(1-p11hat-p12hat-p21hat)/(p12hat*p21hat)
    lnOR <- log(OddsRatio)
    # asymptotic variance of odds ratio Theta
    AsyVarOR <- (1/p11hat + 1/p12hat + 1/p21hat + 1/(1 - p11hat - p12hat - p21hat))
    # test statistics
    z1<-lnOR/sqrt(AsyVarOR/100)
    z2<- Rho/ sqrt(AsyVarRho/100)
    count_lnOR<-ifelse(abs(z1)>1.96,1,0)
    count_Rho<-ifelse(abs(z2)>1.96,1,0)
    # outputs in dataframe by columns
    df <- c(p11hat, p12hat, p21hat, Rho, AsyVarRho, OddsRatio,lnOR, AsyVarOR,count_lnOR,count_Rho)
    cell_counts_data=rbind(df,cell_counts_data)

  }# i

  colnames(cell_counts_data)<-c("p11hat", "p12hat", "p21hat", "Rho", "AsyVarRho", "OddsRatio","lnOR","AsyVarOR","count_lnOR",
                     "count_Rho")
  return(cell_counts_data)
}


result<-replicate(100, expr = {
  Y <- power_function(1000)
  C_lnOR <- as.matrix(table(Y$count_lnOR))
  C_Rho <- as.matrix(table(Y$count_Rho))
  Percent_lnOR <- (C_lnOR[2] / sum(C_lnOR[1], C_lnOR[2])) * 100
  Percent_Rho <- (C_Rho[2] / sum(C_Rho[1], C_Rho[2])) * 100
  c(Percent_lnOR, Percent_Rho)
})

row.names(result)<- c("Power_odds_ratio", "Power_Canonical")
result_transpose<-t(result)

plot(MB,PowerRho, type= "l", col= "red", main= " Power functions " ,
     ylab= "Power", xlab= " Lambda" )
lines(MB,PowerlnTheta,col= "blue" )
legend("bottomright", legend= c("Power_Rho", "Power_lntheta"), col=c("red", "blue"), pch= c(16,17),)

plot(result_transpose$Power_odds_ratio)


