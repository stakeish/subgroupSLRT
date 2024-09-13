rm(list=ls())

library(subgroupSLRT)
library(parallel)
library(abind)
library(nloptr)
library(Rcpp)
library(RcppArmadillo)
load("summary_null.RData")
#source("R/func.R")
#Rcpp::sourceCpp("src/cppshrinkage.cpp")


J <- 3
K <- 3
c <- c(0.05, 2, 2)
B <- 500

nreps <- 2000
gamma.dim <- 10
sample_size <- c(100, 250, 500, 750, 1000)
adj.crt.vec <- summary[, 3]

alpha <- c(1, 2)
beta <- 1
lambda <- 1
gamma <-  rep(1, gamma.dim)
sigma <- 1
alpha.dim <- length(alpha)

summary <- matrix(double(length(sample_size)*3), nrow = length(sample_size), ncol = 3)

for(i in 1:length(sample_size)) {

  n <- sample_size[i]
  adj.crt <- adj.crt.vec[i]

  set.seed(20)

  intercept_set <- array(rep(1, n*nreps), dim = c(n, 1, nreps))
  x_set <- array(rnorm(n*nreps*(alpha.dim - 1)), dim = c(n, alpha.dim - 1, nreps))
  x_set <- abind(intercept_set, x_set, along = 2)
  d_set <- array(rbinom(n*nreps, 1, 0.5), dim = c(n, 1, nreps))
  z_set <- array(rnorm(n*nreps*(gamma.dim - 1)), dim = c(n, gamma.dim - 1, nreps))
  z_set <- abind(intercept_set, z_set, along = 2)
  y_set <- lapply(1:nreps, function(j, x_set, d_set, z_set, alpha, beta, lambda, gamma, sigma){
    x <- x_set[, , j]
    d <- d_set[, , j]
    z <- z_set[, , j]
    n <- nrow(x)
    rlognormal(n, x, d, z, alpha, beta, lambda, gamma, sigma)
  }, x_set = x_set, d_set = d_set, z_set = z_set, alpha = alpha, beta = beta, lambda = lambda, gamma = gamma, sigma = sigma)

  out <- parSapply(cl, 1:nreps, function(j, x_set, d_set, z_set, y_set, J, K, c_bound, B) {

    x <- x_set[, , j]
    d <- d_set[, , j]
    z <- z_set[, , j]
    y <- y_set[[j]]

    c <- c_bound

    out <- bootstrap_test(x, d, z, y, J, c, K, B)

    em.stat <- out$em.stat
    p.val <- out$p.val

    c(em.stat, p.val)
  }, x_set = x_set, d_set = d_set, z_set = z_set, y_set = y_set, J = J, K = K, c_bound = c, B = B)

  em.stat <- out[1, ]
  p.val <- out[2, ]

  rfreq <- mean(p.val < 0.05)
  rfreq.adj <- mean(em.stat > adj.crt)

  summary[i, 1] <- n
  summary[i, 2] <- rfreq
  summary[i, 3] <- rfreq.adj

  save(summary, file = "summary.RData")

  print(n)
}

summary

stopCluster(cl)
