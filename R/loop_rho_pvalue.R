# loop_rho_pvalue
# save.image("one_snp")
# load("one_snp")
# new_data for loop
x <- new_data$x
# loop for all snips
rho <- c(rep(0,111))
p_val<- c(rep(0,111))
for(i in seq(2,112,1)){
  y <- new_data[ ,i]
  rho[i-1] = rho_twoBythree(x,y)
  p_val[i-1] <- exp(-1397*rho[i-1]/2)
}

chi <- c(rep(0,111))
for(i in seq(2,112,1)){
y <- new_data[ ,i]
w = chisq.test(x,y)
chi[i-1] <- w$p.value
}
chi


rho_pval <- data.frame(rho, p_val, chi)
row.names(rho_pval) <- colnames(new_data[-1])
write.table(rho_pval, "rho_pval",
            quote = F,
            row.names = T,
            col.names = T)
which(rho_pval$p_val <= 0.05)
which(chi <= 0.05)
