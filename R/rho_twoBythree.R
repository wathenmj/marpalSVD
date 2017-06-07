rho_twoBythree <- function(x,y){
  temp <- as.matrix(table(x,y))
  tbl <-as.matrix(prop.table(temp))
  col_num <- dim(tbl)[2]
  if(col_num == 1){
    rho <- 0
    return(rho)
  }
  else if(col_num == 2){
    p11hat <-tbl[1,1]; p12hat <- tbl[1,2]
    p21hat <-tbl[2,1]; p22hat <- tbl[2,2]

    p1p <- sum(tbl[1, ]); p2p <- sum(tbl[2, ])
    pp1 <- sum(tbl[, 1]); pp2 <- sum(tbl[,2 ])
    rho <- p11hat^2/(p1p*pp1) + p12hat^2/(p1p*pp2) +
      p21hat^2/(p2p*pp1) + p22hat^2/(p2p*pp2) -
      1
    return(rho)
  }
  else{
    p11hat <-tbl[1,1]; p12hat <- tbl[1,2]; p13hat <- tbl[1,3]
    p21hat <-tbl[2,1]; p22hat <- tbl[2,2]; p23hat <- tbl[2,3]
    p1p <- sum(tbl[1, ]); p2p <- sum(tbl[2, ]) # row marginals
    pp1 <- sum(tbl[, 1]); pp2 <- sum(tbl[,2 ]) # col marginals
    pp3 <- sum(tbl[ ,3]) # col marginals
    b11 <- p11hat/sqrt(p1p*pp1); b12 <- p12hat/sqrt(p1p*pp2); b13 <- p13hat/sqrt(p1p*pp3)
    b21 <- p21hat/sqrt(p2p*pp1); b22 <- p22hat/sqrt(p2p*pp2); b23 <- p23hat/sqrt(p2p*pp3)
    bcells <- c(b11,b12,b13,b21,b22,b23)
    if(sum(bcells)=="NaN"){
      rho <- 0
    }
    else{
      rho <- p11hat^2/(p1p*pp1) + p12hat^2/(p1p*pp2) + p13hat^2/(p1p*pp3) +
        p21hat^2/(p2p*pp1) + p22hat^2/(p2p*pp2) + p23hat^2/(p2p*pp3) - 1
    }
    return(rho)
  }
}


# Quantile 
quantile(DRM.CH, p= seq(0,1,0.1),na.rm=TRUE)
# remove unnecessary column
fms_1=fms[,- c(1,3,53,227:235,237:347)] 
attach(fms_1)
#if the change is less or equal to 0 then 0 otherwise 1
x<-ifelse(DRM.CH <= 0.0,0,1)
table(DRM.CH)
#new data frame
new_data<- data.frame(x,fms_1[ ,1:223] )
dat_y<-new_data[ ,c(2:223)]
dat_x<-new_data[ ,c(1)]

# loop for all snips
rho <- c(rep(0,222))
for(i in 1:222){
  y <- dat_y[,i]
  rho[[i]] = rho_twoBythree(x,y)
}
rho
# Histogram of the rho values
hist(rho, 
     main="Rho 2*3", 
     xlab="Rho", 
     border="blue", 
     col="green",
     xlim=c(-1,0.4),
     las=1, 
     breaks=5)
  

#snips name with rho value
snp<-names(dat_y)
out<-data.frame(snp,rho) 
  



