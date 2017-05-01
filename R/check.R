a <- matrix(sample(0:2, size = 4*3, TRUE), 4, 3) #
a <- prop.table(a)
colnames(a) <- c("geno_1","geno_2","geno_3"); rownames(a) <- c("pop_1","pop_2","pop_3", "pop_4")
ftable(a)
ftable(addmargins(a))

