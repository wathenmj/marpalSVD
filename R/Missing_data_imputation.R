# clean up the data
# multiple imputation 
# 18.2 for NDRM is the cutoff, less than 18.2 is 1, othes 0
#

library(mice)
library(dplyr)
library(ggplot2)
library(genetics)
# Data: we have 223 snips. The NDRM_DIFF are the snips.
#we dont need ID(column=1), ace_id(copumn=3), b2b(column=53), column=227-234,
#column=236-347   

fmsURL<-"http://www.stat-gen.org/book.e1/data/FMS_data.txt"
fms<-read.delim(file=fmsURL, header=TRUE, sep="\t")
attach(fms)
fms_1=fms[,- c(1,3,53,227:235,237:347)]
# subset cases and controls

fms_subset_1 <- fms_1[, -c(28,58:62,97:99, 108, 110, 111, 162:164,186,195)]
dim(fms_subset_1)
#change the level of of different column
levels(fms_1[,28]) <- c("GG","GA","AA","NA")

# Quantile 
quantile(NDRM.CH, p= seq(0,1,0.1),na.rm=TRUE)

x<-ifelse(NDRM.CH <= 18.2,1,0)
cases <-filter(fms_1, x==1)
controls <- filter(fms_1, x == 0)



table(fms_1[,28])
fms_1[,28] <- ifelse(fms_1[,28]== "NA","N/A",fms_1[,28])
which(fms_1[,28]=="Und")
fms_1[135,28]<- NA
fms_1[1391,28]<- NA

fms_1[,58:62]
# column 58 to 62 has only one value level 
table (fms_1[,58])
table (fms_1[,59])
table (fms_1[,60])
table (fms_1[,61])

table (fms_1[,108])
fms_1[,108] <- c("AA","AT","TT","NA")
table (fms_1[,108])

# please check this 
table (fms_1[,110])

table (fms_1[,111])

# need o remove this
table (fms_1[,162])
#check this 
table (fms_1[,163])


# work with cases first
attach(cases)
cases[1:6,1:5]


###miss percentage by row and column
pMiss <- function(x){sum(is.na(x))/length(x)*100}
apply(cases,2,pMiss)
apply(cases,1,pMiss)
apply(controls,2,pMiss)

#see the missing pattern
library(VIM)
aggr_plot <- aggr(cases, col=c('navyblue','red'), numbers=TRUE, 
                  sortVars=TRUE, labels=names(cases), cex.axis=.7, gap=3,
                  ylab=c("Histogram of missing data","Pattern"))


aggr_plot <- aggr(controls, col=c('navyblue','red'), numbers=TRUE, sortVars=TRUE, 
                  labels=names(controls), cex.axis=.7, gap=3, 
                  ylab=c("Histogram of missing data","Pattern"))

#imputing missing data
tempDaat <- mice(cases,m=1,seed=500) # this is case 
summary(tempData)

tempDaat1<-mice(controls,m=1,seed=500) # this is control 

#check the missing data imputation by variable 
tempData$imp$varable_name # variable_name will be snip name

# to check the imputation method
tempDaat$meth

# to check the complete cases again 
completedData <- complete(tempData,1)

# to check the distribution of the 
densityplot(tempData)
