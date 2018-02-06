# chi_square_against_Rho

indep <- c(0.06,0.3,0.24,0.04,0.2,0.16)
t_1 <- matrix(indep,nrow = 2, ncol = 3, byrow = TRUE)


a <- t_1[1,1]; b <- t_1[1,2]; c <- t_1[1,3]
d <- t_1[2,1]; e <- t_1[2,2]; f <- t_1[2,3]

x <- rchisq(1000, 5/2)

hist(t_2, prob=TRUE, xlab = "")
curve(dchisq(x, df=2), col='green', add=TRUE)
curve(dchisq(x, df=5/2), col='red', add=TRUE )
mtext("Red df = 5/2 Green df = 2", side = 3)
mtext("Distribution of 1000 Rhos for 100 samples", side = 1, line = 2.5)
mtext("Sample <- rmultinom(1000, 100, prob= c(0.06,0.3,0.24,0.04,0.2,0.16))", side = 1, line =4)




chisq.test(t_1,5/2)
t_1
chisq.test(t_1,5/2)
