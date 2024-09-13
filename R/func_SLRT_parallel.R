#' functions required to implement SLRT

#' A function to generate a set of initial values for beta, lambda and gamma
#' @param ninits the number of initial values
#' @gamma.dim the dimension of gamma
#' @alpha.dim the dimension of alpha
#' @gamma.bound bound for gamma


shrinkageMLEinit <- function(ninits, gamma.dim, alpha.dim, gamma.bound = c(-2, 2)){
  alpha <- matrix(runif(alpha.dim*ninits, min = -3, max = 3), nrow=alpha.dim)
  beta1 <- matrix(runif(ninits, min = -3, max=3), nrow = 1)
  beta2 <- matrix(runif(ninits, min = -3, max=3), nrow= 1)
  sigma <- matrix(runif(ninits, min = 0.5, max=5), nrow=1)
  gamma.lower <- gamma.bound[1]
  gamma.upper <- gamma.bound[2]
  gamma <- matrix(runif(ninits*gamma.dim, min=gamma.lower, max=gamma.upper), nrow=gamma.dim)

  list(alpha = alpha, beta1=beta1, beta2=beta2, gamma=gamma, sigma = sigma)
}

#' A funciton to calculate shrinkage likelihood ratio test statistics
#' @param x the regressor matrix X
#' @param d the treatment indicator vector D
#' @param z the regressor matrix for classification Z
#' @param y the outcome vector y
#' @param p the value of tuning parameter p
#' @param ninits the number of initial values for estimation of smle
#' @param ninits.equal the number of initial values for estimation of mle with equal weight
#' @param epsilon threshold value for EM algorithm
#' @param maxit the maximum number of iterations for EM algorithm
#' @param epsilon.short threshold value for short EM algorithm
#' @param maxit.short the maximum number of iterations for short EM algorithm

SLRT_parallel <- function(x, d, z, y, p, ninits, ninits.equal, epsilon = 1e-08, maxit = 500, epsilon.short = 1e-02, maxit.short = 100, n.core = detectCores() - 1) {

  set.seed(123456)
  q1 <- ncol(x)
  n <- length(y)
  xx <- t(x) %*% x
  xy <- t(x) %*% y

  # Model estimation under the null
  out.null <- homogeneousMLE(x, d, y)
  ll.null <- out.null$ll

  # Model estimation under the alternative with gamma = 0 (pi = 0.5)
  y_0 <- y[d == 0]
  x_0 <- x[d == 0, ]
  alpha.0 <- solve(t(x_0) %*% x_0, t(x_0) %*% y_0)
  sigma.0 <- sqrt(mean((y_0 - x_0 %*% alpha.0)^2))

  alpha.initial <- matrix(rep(alpha.0, ninits.equal), nrow = q1, ncol = ninits.equal)
  beta1.initial <- runif(ninits.equal, min = -5, max = 5)
  beta2.initial <- runif(ninits.equal, min = -5, max = 5)
  sigma.init <- rep(sigma.0, ninits.equal)

  out.equal.cpp <- cppem_normal(y, x, z, d,
                    0.5, alpha.initial, beta1.initial, beta2.initial,
                    sigma.0, ninits.equal, maxit, epsilon)

  alpha.set <- out.equal.cpp$alpha
  beta1.set <- c(out.equal.cpp$beta1)
  beta2.set <- c(out.equal.cpp$beta2)
  sigma.set <- c(out.equal.cpp$sigma)
  ll.set <- c(out.equal.cpp$ll)
  rm(out.equal.cpp)

  index.equal <- which.max(ll.set)
  alpha.equal <- alpha.set[, index.equal]
  beta1.equal <- beta1.set[index.equal]
  beta2.equal <- beta2.set[index.equal]
  sigma.equal <- sigma.set[index.equal]
  ll.equal <- ll.set[index.equal]
  rm(ll.set, alpha.set, beta1.set, beta2.set, sigma.set)

  # Maximum shrinkage likelihood estimation
  q2 <- ncol(z)
  ninits.short <- ninits*q2

  tmp <- shrinkageMLEinit(ninits = ninits.short, gamma.dim = q2, alpha.dim = q1)

  theta0 <- rbind(tmp$alpha, tmp$beta1, tmp$beta2, tmp$sigma, tmp$gamma)
  theta0[, 1] <- c(alpha.equal, beta1.equal, beta2.equal, sigma.equal, double(q2))
  theta0.list <- lapply(seq_len(ncol(theta0)), function(i) theta0[,i])
  rm(tmp, theta0)

  out.short <- mclapply(theta0.list, function(theta.init, x, d, z, y, p, epsilon, maxit, xx, xy, q1, q2){

    out.em <- em_shrinkage(x, d, z, y, p, epsilon, maxit, theta.init, xx, xy, q1, q2)

    out.em

  }, x = x, d = d, z = z, y = y, p = p, epsilon = epsilon.short, maxit = maxit.short, xx = xx, xy = xy, q1 = q1, q2 = q2, mc.cores = n.core)

  ploglik.set <- sapply(out.short, function(x) x$penloglik)
  alpha.set <- sapply(out.short, function(x) x$alpha)
  beta1.set <- sapply(out.short, function(x) x$beta1)
  beta2.set <- sapply(out.short, function(x) x$beta2)
  sigma.set <- sapply(out.short, function(x) x$sigma)
  gamma.set <- sapply(out.short, function(x) x$gamma)
  theta.short.set <- rbind(alpha.set, beta1.set, beta2.set, sigma.set, gamma.set)
  rm(out.short, alpha.set, beta1.set, beta2.set, sigma.set, gamma.set)

  index.ninits <- order(ploglik.set, decreasing = TRUE)[1:ninits]
  theta.short.set <- theta.short.set[, index.ninits]
  theta.short.list <- lapply(seq_len(ncol(theta.short.set)), function(i) theta.short.set[,i])
  rm(theta.short.set)


  out.long <- mclapply(theta.short.list, function(theta.init, x, d, z, y, p, epsilon, maxit, xx, xy, q1, q2){

    out.em <- em_shrinkage(x, d, z, y, p, epsilon, maxit, theta.init, xx, xy, q1, q2)

    out.em

  }, x = x, d = d, z = z, y = y, p = p, epsilon = epsilon, maxit = maxit, xx = xx, xy = xy, q1 = q1, q2 = q2, mc.cores = n.core)

  ploglik.set <- sapply(out.long, function(x) x$penloglik)
  index.top <- which.max(ploglik.set)

  alpha.hat <- out.long[[index.top]]$alpha
  beta1.hat <- out.long[[index.top]]$beta1
  beta2.hat <- out.long[[index.top]]$beta2
  sigma.hat <- out.long[[index.top]]$sigma
  gamma.hat <- out.long[[index.top]]$gamma
  ploglik <- out.long[[index.top]]$penloglik
  rm(out.long)
  loglik <- ploglik + p*sum(abs(gamma.hat))

  par <- list(alpha=alpha.hat, beta1=beta1.hat, beta2=beta2.hat, sigma=sigma.hat, gamma=gamma.hat)

  SLRT <- 2*(max(ll.equal, loglik) - ll.null)
  LRT.equal <- 2*(ll.equal - ll.null)

  list(SLRT = SLRT, LRT.equal = LRT.equal, par = par)

}

#' A funciton to implement EM-algorithm
#' @param x the regressor matrix X
#' @param d the treatment indicator vector D
#' @param z the regressor matrix for classification Z
#' @param y the outcome vector y
#' @param p the value of tuning parameter p
#' @param epsilon threshold value for EM algorithm
#' @param maxit the maximum number of iterations for EM algorithm
em_shrinkage <- function(x, d, z, y, p, epsilon, maxit, theta.init, xx, xy, q1, q2){

  alpha.k <- theta.init[1:q1]
  beta1.k <- theta.init[q1+1]
  beta2.k <- theta.init[q1+2]
  sigma.k <- theta.init[q1+3]
  gamma.k <- theta.init[(q1+4):(q1+q2+3)]
  diff <- 1.0
  oldpenloglik <- -Inf
  x_alpha <- x %*% alpha.k
  z_gamma <- z %*% gamma.k
  w.mat <- matrix(nrow=n, ncol=2)
  #xx <- t(x) %*% x
  #xy <- t(x) %*% y
  mean.1 <- y - x_alpha - d*beta1.k
  mean.2 <- y - x_alpha - d*beta2.k

  for(jj in 1:maxit){

    ##E-step
    #mean.1 <- y - x_alpha - d*beta1.k
    #mean.2 <- y - x_alpha - d*beta2.k
    log.prop <- plogis(z_gamma)
    w.mat[, 1] <- log.prop * dnorm(mean.1, 0, sigma.k)
    w.mat[, 2] <- (1-log.prop) * dnorm(mean.2, 0, sigma.k)
    l.vec <- rowSums(w.mat)
    w.mat <- w.mat / l.vec

    ##Exit loop or not
    penloglik <- sum(log(l.vec)) - p*norm(as.matrix(gamma.k))
    diff <- penloglik - oldpenloglik
    oldpenloglik <- penloglik

    #print(diff)

    #if(diff < 0 || penloglik == -Inf) {
    #  if(penloglik == -Inf){
    #    print("penloglik is -Inf")
    #  } else {
    #    print("diff < 0")
    #  }
    #  break
    #}

    if(diff < epsilon || penloglik == -Inf){
      break
    }

    ##M-step
    w1 <- w.mat[, 1]
    w2 <- w.mat[, 2]
    w1_d <- w1 * d
    w2_d <- w2 * d
    x_w1_d <- t(x) %*% w1_d
    x_w2_d <- t(x) %*% w2_d
    mat.equation.alpha <- cbind(xx, x_w1_d, x_w2_d)
    mat.equation.beta1 <- cbind(t(x_w1_d), sum(w1_d), 0)
    mat.equation.beta2 <- cbind(t(x_w2_d), 0, sum(w2_d))
    mat.equation <- rbind(mat.equation.alpha, mat.equation.beta1, mat.equation.beta2)

    vec.equation.alpha <- xy
    vec.equation.beta1 <- sum(y * w1_d)
    vec.equation.beta2 <- sum(y * w2_d)
    vec.equation <- c(vec.equation.alpha, vec.equation.beta1, vec.equation.beta2)

    coef.k <- try(solve(mat.equation, vec.equation), silent = TRUE)
    if(class(coef.k) == "try-error"){
      break
    }

    alpha.k <- coef.k[1:q1]
    beta1.k <- coef.k[q1+1]
    beta2.k <- coef.k[q1+2]
    x_alpha <- x %*% alpha.k

    mean.1 <- y - x_alpha - d*beta1.k
    mean.2 <- y - x_alpha - d*beta2.k

    sigma.k <- sqrt(mean(w1 * (mean.1)^2) + mean(w2 * (mean.2)^2))

    out.logit <- lbfgs::lbfgs(llogit, dlogit, gamma.k, invisible=1, orthantwise_c = p, z=z, w=w1)
    gamma.k <- out.logit$par
    z_gamma <- z %*% gamma.k

  }

  list(alpha = alpha.k, beta1 = beta1.k, beta2 = beta2.k, sigma = sigma.k, gamma = gamma.k, penloglik = penloglik)

}
