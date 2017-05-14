# url_non_dominant_arm_data
fmsURL<-"http://www.stat-gen.org/book.e1/data/FMS_data.txt"
fms<-read.delim(file=fmsURL, header=TRUE, sep="\t")
attach(fms)

write.table(fms,"fms",
            quote = F,
            row.names = F,
            col.names = T)

# see page 21 of Applied Statistical Genetics with R. Andrea S. Foulkes
colnames(fms)
GenoCount <- summary(actn3_rs540874)
GenoCount

NumbObs <- sum(!is.na(actn3_rs540874))
NumbObs

# genotype frequencies for AA, GA, GG, and NA's are given respectively
GenoFreq <- as.vector(GenoCount/NumbObs)
GenoFreq

# frequencies  of A and G alleles are calulated as follows

FreqA <- (2*GenoFreq[1] + GenoFreq[2])/2
FreqA

FreqG <- (2*GenoFreq[3] + GenoFreq[2])/2
FreqG
# so A is the minor Allele with a frequency of 0.431

library(genetics); library(coin)
# Cochran-Armitage (C-A) trend test p.42
Geno <- genotype(actn3_rs540874, sep = "")
summary(Geno)

Geno <- esr1_rs1042717
Trait <- as.numeric(pre.BMI>25)

GenoOrd <- ordered(Geno)
independence_test(Trait~GenoOrd, teststat ="quad",
                  scores=list(GenoOrd=c(0,1,2)))
