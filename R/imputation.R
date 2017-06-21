# imputation
setwd("C:/marpalSVD/R")
# load("SVD")
# save.image("SVD")
library(mice)

# Following is in the image space; no need to run ####
Data_imput <- mice(fms_subset_2,m=1,seed=500)
Data_imput_data <- complete(Data_imput,action=1,include=TRUE)
x<-ifelse(Data_imput_data$NDRM.CH <= 18.2,1,0)
x<-data.frame(Data_imput_data$NDRM.CH)
y<-data.frame(Data_imput_data [,1:111])
new_data<- data.frame(x,y )
dat_y<-new_data[ ,c(2:112)]
dat_x<-new_data[ ,c(1)]

# loop for all snips ####
rho <- c(rep(0,111))
for(i in 1:111){
  x<-dat_x
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
     xlim=c(-0.5,0.5),
     las=1,
     breaks=1)
