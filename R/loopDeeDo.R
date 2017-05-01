snpNo <- dim(geno_example)[2] - 2  # first two columns are ids and pops
RhoChr <- c(rep(0,snpNo))
timestamp() #takes about 45 minutes to run, try subsetting to get some faster results
for (i in 1:snpNo) {
  a <- geno_example[ ,2]
  b <- geno_example[,i+2]
  Rho(a,b) -> RhoChr[i]
}
RhoChr
