rm(list=ls())

library(subgroupSLRT)
library(abind)
library(lbfgs)
library(parallel)

load("coefficients.RData")
intercept <- coefficients[1]
coef <- coefficients[2]

gamma.dim <- 100
n <- 1000

alpha <- c(1, 2)
beta <- 1
lambda <- 1
gamma <-  rep(1, gamma.dim)
sigma <- 1
alpha.dim <- length(alpha)

p <- intercept + coef * n^(7/8) * sqrt(log(gamma.dim))

set.seed(20)
intercept <- rep(1, n)
x <- matrix(rnorm(n*alpha.dim), nrow=n)
d <- rbinom(n, 1, 0.5)
z <- matrix(rnorm(n*(gamma.dim-1)), nrow=n)
z <- cbind(intercept, z)
y <- rlognormal(n, x, d, z, alpha, beta, lambda, gamma, sigma)

time_start <- proc.time()
out <- SLRT_parallel(x, d, z, y, p, ninits = 10, ninits.equal = 10)
time_end <- proc.time()
time_elapsed <- time_end[3] - time_start[3]
sink("time_elapsed.txt")
print(time_elapsed)
sink()
save(time_elapsed, file="time_elapsed.RData")
