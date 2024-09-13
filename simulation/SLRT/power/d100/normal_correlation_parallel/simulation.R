rm(list=ls())
gc();  gc();

library(subgroupSLRT)
library(parallel)
library(abind)
library(lbfgs)

load("coefficients.RData")
load("summary_null.RData")
load("LRT.equal_vec.RData")
load("SLRT_vec.RData")

intercept <- coefficients[1]
coef <- coefficients[2]
nreps <- 5000
gamma.dim <- 100
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


set.seed(1234)
Sigma_mat_root <- matrix(runif((gamma.dim - 1)^2, min = -1, max = 1), ncol = gamma.dim - 1)
Sigma_mat_square <- Sigma_mat_root^2
Sigma_mat_sum <- apply(Sigma_mat_square, 1, sum)
Sigma_mat_root <- Sigma_mat_root / sqrt(Sigma_mat_sum)
Sigma_mat <- Sigma_mat_root %*% t(Sigma_mat_root)


for(i in 5:5) {

  n <- sample_size[i]
  p <- intercept + coef * n^(7/8) * sqrt(log(gamma.dim))
  adj_crt_equal <- adj_crt_vec_equal[i]
  adj_crt_SLR <- adj_crt_vec_SLR[i]


  set.seed(20)

  intercept_set <- array(rep(1, n*nreps), dim = c(n, 1, nreps))
  x_set <- array(rnorm(n*nreps*(alpha.dim - 1)), dim = c(n, alpha.dim - 1, nreps))
  x_set <- abind(intercept_set, x_set, along = 2)
  d_set <- array(rbinom(n*nreps, 1, 0.5), dim = c(n, 1, nreps))
  z_set <- array(unlist(lapply(1:nreps, function(j, Sigma_mat, gamma.dim, n){
    norm_part <- MASS::mvrnorm(n, mu = rep(0, gamma.dim - 1), Sigma = Sigma_mat)
    cbind(rep(1, n), norm_part)
  }, Sigma_mat = Sigma_mat, gamma.dim = gamma.dim, n = n)), dim = c(n, gamma.dim, nreps))
  #z_set_big <- array(unlist(z_set), dim = c(n, gamma.dim, nreps))
  y_set <- lapply(1:nreps, function(j, x_set, d_set, z_set, alpha, beta, lambda, gamma, sigma){
    x <- x_set[, , j]
    d <- d_set[, , j]
    z <- z_set[, , j]
    n <- nrow(x)
    rlognormal(n, x, d, z, alpha, beta, lambda, gamma, sigma)
  }, x_set = x_set, d_set = d_set, z_set = z_set, alpha = alpha, beta = beta, lambda = lambda, gamma = gamma, sigma = sigma)
  rm(intercept_set)
  gc()

  #SLRT_vec <- double(nreps)
  #LRT.equal_vec <- double(nreps)

  for(iter in 5000:nreps){
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

  rfreq_slrt <- mean(SLRT_vec > critical_value)
  rfreq_slrt_adj <- mean(SLRT_vec > adj_crt_SLR)
  rfreq_lrt.equal <- mean(LRT.equal_vec > critical_value)
  rfreq_lrt.equal_adj <- mean(LRT.equal_vec > adj_crt_equal)

  summary_matrix[i, ] <- c(rfreq_lrt.equal, rfreq_slrt, rfreq_lrt.equal_adj,  rfreq_slrt_adj)

  save(summary_matrix, file = "summary.RData")

  print(n)

}

summary_matrix
