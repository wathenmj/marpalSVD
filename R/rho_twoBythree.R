rho_twoBythree <- function(x,y){
  temp <- as.matrix(table(x,y))
  tbl <-as.matrix(prop.table(temp))
  col_num <- dim(tbl)[2]
  if(col_num == 1){
    rho <- 0
    return(rho)
  }
  else if(col_num == 2){
    p11 <-tbl[1,1]; p12 <- tbl[1,2]
    p21 <-tbl[2,1]; p22 <- tbl[2,2]

    p1p <- sum(tbl[1, ]); p2p <- sum(tbl[2, ])
    pp1 <- sum(tbl[, 1]); pp2 <- sum(tbl[,2 ])
    rho <- p11^2/(p1p*pp1) + p12^2/(p1p*pp2) +
      p21^2/(p2p*pp1) + p22^2/(p2p*pp2) -
      1
    return(rho)
  }
  else{
    p11 <-tbl[1,1]; p12 <- tbl[1,2]; p13 <- tbl[1,3]
    p21 <-tbl[2,1]; p22 <- tbl[2,2]; p23 <- tbl[2,3]
    p1p <- sum(tbl[1, ]); p2p <- sum(tbl[2, ]) # row marginals
    pp1 <- sum(tbl[, 1]); pp2 <- sum(tbl[,2 ]) # col marginals
    pp3 <- sum(tbl[ ,3]) # col marginals
    b11 <- p11/sqrt(p1p*pp1); b12 <- p12/sqrt(p1p*pp2); b13 <- p13/sqrt(p1p*pp3)
    b21 <- p21/sqrt(p2p*pp1); b22 <- p22/sqrt(p2p*pp2); b23 <- p23/sqrt(p2p*pp3)
    bcells <- c(b11,b12,b13,b21,b22,b23)
    if(sum(bcells)=="NaN"){
      rho <- 0
    }
    else{
      rho <- p11^2/(p1p*pp1) + p12^2/(p1p*pp2) + p13^2/(p1p*pp3) +
        p21^2/(p2p*pp1) + p22^2/(p2p*pp2) + p23^2/(p2p*pp3) - 1
    }
    return(rho)
  }
}






