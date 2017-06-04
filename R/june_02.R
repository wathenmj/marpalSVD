# june_02
#fmsURL<-"http://www.stat-gen.org/book.e1/data/FMS_data.txt"
#fms<-read.delim(file=fmsURL, header=TRUE, sep="\t")
#attach(fms)
# fms1 <- fms[ , -c(1, 3, 53, 227:234, 236:347)]
# dim(fms1) # dim 1397 347
# write.table(fms1, "fms_subset_1",
  #          quote = F,
   #         row.names = F,
    #        col.names = T)
fms_subset_1 <- read.delim("C:/marpalSVD/fms_subset_1", sep = "")
attach(fms1)
x <- ifelse(NDRM_DIFF <= 3.8, 1, 0) #make case at tenth percentile
# x will be first argument of rho_twoBythree
table(x)
# load function rho_twoBythree.R
# test rho_twoBythree with example
fms_subset_1[1:2,1:6]
y <- actn3_r577x
rho_twoBythree(x,y) # gives output value [1] 0.01187055
