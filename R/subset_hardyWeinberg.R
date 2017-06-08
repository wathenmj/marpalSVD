# subset_hardyWeinberg
library(dplyr);library(genetics)
fms_subset_1 <- read.delim("fms_subset_1", sep = "")
attach(fms_subset_1)
x <- ifelse(NDRM_DIFF <= 3.8, 1, 0) #make case at tenth percentile
# x will be first argument of rho_twoBythree

# subset cases and controls
cases <- filter(fms_subset_1, NDRM_DIFF <= 3.8)
controls <- filter(fms_subset_1, NDRM_DIFF > 3.8)

# work with cases first
attach(cases)
cases[1:6,1:5]

snp1 <- actn3_rs540874 #testing one example snp
ObsCount <- table(snp1);ObsCount
Nobs <- sum(ObsCount);Nobs
FreqG <- (2*ObsCount[3] + ObsCount[2])/(2*Nobs);FreqG
ExpCount <- c(Nobs*(1 - FreqG^2), 2*Nobs*FreqG*(1-FreqG), Nobs*FreqG^2);ExpCount
ChiSqStat <- sum((ObsCount - ExpCount)^2/ExpCount);ChiSqStat
qchisq(1-0.5, df = 1)
# Since 35.49683 > 0.4549364 do not reject the null; the alleles on the two homologous
#   chromosomes are not associated with each other

# try the same example with controls instead
attach(controls)
snp1 <- actn3_rs540874 #testing one example snp
ObsCount <- table(snp1);ObsCount
Nobs <- sum(ObsCount);Nobs
FreqG <- (2*ObsCount[3] + ObsCount[2])/(2*Nobs);FreqG
ExpCount <- c(Nobs*(1 - FreqG^2), 2*Nobs*FreqG*(1-FreqG), Nobs*FreqG^2);ExpCount
ChiSqStat <- sum((ObsCount - ExpCount)^2/ExpCount);ChiSqStat
qchisq(1-0.5, df = 1)

# Since 315.0774 > 0.4549364 do not reject the null; the alleles on the two homologous
#   chromosomes are not associated with each other
# Use genetics R package
attach(cases) # work with cases only
snp1 <- genotype(actn3_rs540874, sep = "")
HWE.chisq(snp1)
# again fail to reject; mating appears to be random for cases
attach(controls) # work with cases only
snp1 <- genotype(actn3_rs540874, sep = "")
HWE.chisq(snp1)
# random mating appears to hold for controls as well
# genetics R package has Fishers Exact test
#   for use when expected cell counts fall below below 5
HWE.exact(snp1)

snp <- rep(0,207)
hwe_pvalue <- rep(0,207)
expect_count <- rep(0,207)
exact_pvalue <- rep(0,207)

fms_subset_1[,28]# Levels: AA GA GG Und
fms_subset_1[,58]# Levels: GG
fms_subset_1[ ,59] #Levels: CC
for(i in 1:207){
  snp <- genotype(fms_subset_1[,i], sep = "")
  k <- HWE.chisq(snp)
  k$p.value -> hwe_pvalue[i]
  l <- HWE.exact(snp)
  l$p.value -> exact_pvalue[i]
}

fms_subset_1 <- fms_subset_1[, -c(28,58:62,97:99, 108, 110, 111, 162:164,186,195)]
