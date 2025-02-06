rm(list=ls())

library(subgroupSLRT)
library(abind)
library(lbfgs)
library(parallel)

gamma.dim <- 100
n <- 1000

J <- 160
K <- 6
c <- c(0.05, 2, 2)
B <- 500
which.parallel <- "J"

alpha <- c(1, 2)
beta <- 1
lambda <- 1
gamma <-  rep(1, gamma.dim)
sigma <- 1
alpha.dim <- length(alpha)


set.seed(20)
intercept <- rep(1, n)
x <- matrix(rnorm(n*alpha.dim), nrow=n)
d <- rbinom(n, 1, 0.5)
z <- matrix(rnorm(n*(gamma.dim-1)), nrow=n)
z <- cbind(intercept, z)
y <- rlognormal(n, x, d, z, alpha, beta, lambda, gamma, sigma)

time_start <- proc.time()
out <- bootstrap_test_wo_parallel(x=x, d=d, z=z, y=y, J=J, c=c, K=K, B=B, which.parallel = which.parallel)
time_end <- proc.time()
time_elapsed <- time_end[3] - time_start[3]
sink("time_elapsed.txt")
print(time_elapsed)
sink()
save(time_elapsed, file="time_elapsed.RData")
