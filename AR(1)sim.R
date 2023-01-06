Replication <- 1e3
N.Perm <- 1e4
Proportion <- rep(NA, Replication)
for(reps in 1 : Replication){
#generate the data
set.seed(23 *reps)
treatment <- arima.sim(list(order=c(1,0,0), ar=.1), n=1e3)
control <- arima.sim(list(order=c(1,0,0), ar=.1), n=1e3)
#Run the permutation test
Perm <- mc_permutation_test(treatment = treatment, control = control, nsim = N.Perm)
## Forth, decide whether to reject or not
Proportion[reps] <- mean(abs(Perm$perm) > abs(Perm$obs))

## Check simulation progress
cat("Replication =", reps, "\n")
}

hist(Proportion)
mean(Proportion <= 0.05) #How many times the test was rejected
beep(5)



