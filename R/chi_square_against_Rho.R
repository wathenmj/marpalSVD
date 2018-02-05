indep <- c(0.06,0.3,0.24,0.04,0.2,0.16)
t_1 <- matrix(indep,nrow = 2, ncol = 3, byrow = TRUE)
t_1 <- prop.table(t_1)

a <- t_1[1,1]; b <- t_1[1,2]; c <- t_1[1,3]
d <- t_1[2,1]; e <- t_1[2,2]; f <- t_1[2,3]

x <- rchisq(1000, 2)

hist(t_2, prob=TRUE)
curve(dchisq(x, df=2), col='green', add=TRUE)
curve(dchisq(x, df=5/2), col='red', add=TRUE )

chisq.test(t_1,5/2)
t_1
chisq.test(t_1,5/2)
