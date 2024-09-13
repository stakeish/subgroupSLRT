rm(list=ls())

library(subgroupSLRT)
library(parallel)
library(abind)
library(nloptr)
library(Rcpp)
library(RcppArmadillo)
load("summary_null.RData")

J <- 10
K <- 6
c <- c(0.05, 2, 2)
B <- 500

nreps <- 2000
gamma.dim <- 100
sample_size <- c(100, 250, 500, 750, 1000)
adj.crt.vec <- summary_matrix[, 3]

alpha <- c(1, 2)
beta <- 1
lambda <- 1
gamma <-  rep(1, gamma.dim)
sigma <- 1
alpha.dim <- length(alpha)

summary_matrix <- matrix(double(length(sample_size)*3), nrow = length(sample_size), ncol = 3)

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
  rm(intercept_set)

  em_vec <- double(nreps)
  pval_vec <- double(nreps)

  for(iter in 1:nreps){
    x <- x_set[, , iter]
    d <- d_set[, , iter]
    z <- z_set[, , iter]
    y <- y_set[[iter]]
    print(iter)
    out <- bootstrap_test(x, d, z, y, J, c, K, B, n.init = 10)
    em_vec[iter] <- out$em.stat
    pval_vec[iter] <- out$p.val
    if(iter %% 100 == 0) {
      print(c(n, iter))
    }
  }

  rfreq <- mean(pval_vec < 0.05)
  rfreq.adj <- mean(em_vec > adj.crt)

  summary_matrix[i, 1] <- n
  summary_matrix[i, 2] <- rfreq
  summary_matrix[i, 3] <- rfreq.adj

  save(summary_matrix, file = "summary.RData")

}

summary_matrix
