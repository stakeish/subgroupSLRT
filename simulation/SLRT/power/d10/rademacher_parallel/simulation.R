rm(list=ls())

library(subgroupSLRT)
library(parallel)
library(abind)
library(lbfgs)

load("coefficients.RData")
load("summary_null.RData")
intercept <- coefficients[1]
coef <- coefficients[2]
nreps <- 5000
gamma.dim <- 10
sample_size <- c(100, 250, 500, 750, 1000)

adj_crt_vec_equal <- summary_matrix[, 3]
adj_crt_vec_SLR <- summary_matrix[, 4]
critical_value <- qnorm(0.95)^2

alpha <- c(1, 2)
beta <- 1
lambda <- 1
gamma <-  rep(1, gamma.dim)
sigma <- 1
alpha.dim <- length(alpha)

summary_matrix <- matrix(double(length(sample_size)*4), nrow = length(sample_size), ncol = 4)
colnames(summary_matrix) <- c("power equal", "power SLR", "adj. power equal", "adj. power SLR")

for(i in 1:length(sample_size)) {

  n <- sample_size[i]
  p <- intercept + coef * n^(7/8) * sqrt(log(gamma.dim))
  adj_crt_equal <- adj_crt_vec_equal[i]
  adj_crt_SLR <- adj_crt_vec_SLR[i]

  set.seed(20)

  intercept_set <- array(rep(1, n*nreps), dim = c(n, 1, nreps))
  x_set <- array(rnorm(n*nreps*(alpha.dim - 1)), dim = c(n, alpha.dim - 1, nreps))
  x_set <- abind(intercept_set, x_set, along = 2)
  d_set <- array(rbinom(n*nreps, 1, 0.5), dim = c(n, 1, nreps))
  z_set <- array((-1)^rbinom(n*nreps*(gamma.dim - 1), 1, 0.5), dim = c(n, gamma.dim - 1, nreps))
  z_set <- abind(intercept_set, z_set, along = 2)
  y_set <- lapply(1:nreps, function(j, x_set, d_set, z_set, alpha, beta, lambda, gamma, sigma){
    x <- x_set[, , j]
    d <- d_set[, , j]
    z <- z_set[, , j]
    n <- nrow(x)
    rlognormal(n, x, d, z, alpha, beta, lambda, gamma, sigma)
  }, x_set = x_set, d_set = d_set, z_set = z_set, alpha = alpha, beta = beta, lambda = lambda, gamma = gamma, sigma = sigma)
  rm(intercept_set)

  SLRT_vec <- double(nreps)
  LRT.equal_vec <- double(nreps)

  for(iter in 1:nreps){
    x <- x_set[, , iter]
    d <- d_set[, , iter]
    z <- z_set[, , iter]
    y <- y_set[[iter]]
    print(iter)
    out <- SLRT_parallel(x, d, z, y, p, ninits = 10, ninits.equal = 10)
    SLRT_vec[iter] <- out$SLRT
    LRT.equal_vec[iter] <- out$LRT.equal
    save(SLRT_vec, file = "SLRT_vec.RData")
    save(LRT.equal_vec, file = "LRT.equal_vec.RData")
    if(iter %% 1000 == 0) {
      print(c(n, iter))
    }
  }

  rm(x_set, d_set, z_set, y_set)

  rfreq_slrt <- mean(SLRT_vec > critical_value)
  rfreq_slrt_adj <- mean(SLRT_vec > adj_crt_SLR)
  rfreq_lrt.equal <- mean(LRT.equal_vec > critical_value)
  rfreq_lrt.equal_adj <- mean(LRT.equal_vec > adj_crt_equal)

  summary_matrix[i, ] <- c(rfreq_lrt.equal, rfreq_slrt, rfreq_lrt.equal_adj,  rfreq_slrt_adj)

  save(summary_matrix, file = "summary.RData")

  print(i)

}

summary_matrix

stopCluster(cl)
