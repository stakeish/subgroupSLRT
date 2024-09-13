#' functions required for general purposes

#' A function to generate y from a logistic-normal mixture model
#' @param n sample size
#' @param x design matrix x including intercept
#' @param d vector of treatment indicator
#' @param z design matrix z for subgroup classification
#' @param alpha the value of alpha
#' @param beta the value of beta
#' @param lambda the value of lambda
#' @param gamma the value of gamma
#' @param sigma the value of sigma

rlognormal <- function(n, x, d, z, alpha, beta, lambda, gamma, sigma){
  membership <- rbinom(n, 1, plogis(z %*% gamma))
  subgroup_mean <- ifelse(membership==1, x %*% alpha + d * (beta + lambda), x %*% alpha + d * beta)
  rnorm(n, subgroup_mean, sigma)
}

#' A function to implement MLE under the null model
#' @param y vecotr of outcome variable y

homogeneousMLE <- function(x, d, y){
  x.dim <- ncol(x)
  X <- cbind(x, d)
  estimator <- solve(t(X) %*% X) %*% t(X) %*% y
  alpha.hat <- estimator[1:x.dim]
  beta.hat <- estimator[x.dim+1]
  sigma.hat <- sqrt(mean((y - X %*% estimator)^2))
  ll <- sum(log(dnorm(y - X %*% estimator, mean=0, sd=sigma.hat)))
  par <- list(alpha = alpha.hat, beta = beta.hat, sigma = sigma.hat)

  list(ll = ll, par=par)
}
