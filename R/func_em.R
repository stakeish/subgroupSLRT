#' functions required to conduct EM test proposed by Shen and He (2015) JASA

#' A function to conduct EM test
#' @param x confounding covariates used for the outcome equation
#' @param d treatment variables
#' @param z covariates used for subgroup classification including intercept
#' @param J the number of initial gammas
#' param gamma_space the parameter space for gamma
#' @param c prespecified threshold values for parameter space of gamma
#' @param K the number of EM iteration
#' @param b the number of bootstrap replications
#' @param n.init the number of initial values for EM algorithm for MLE with gamma fixed
#' @param upper upper value for the parameter space of alpha, beta, lambda and sigma
#' @param lower lower value for the parameter space of alpha and beta
#' @param lower.zero lower positive value for the parameter space of lambda and sigma
#' @param epsilon threshold value for EM algorithm
#' @param maxit the maximum number of iterations for EM algorithm

em_test <- function(x, d, z, y, J, c, K, xx, xy, n.init = 20, upper = 5, lower = -5, lower.zero = 0.05, epsilon = 1e-08, maxit = 500, n.core = detectCores() - 1) {

  n <- length(y)
  q2 <- ncol(z)
  gamma.init.set <- generate_gamma(c, q2, J)
  q1 <- ncol(x)

  theta.init.set <- theta_initial(q1, upper, lower, lower.zero, n.init)

  out <- lapply(1:J, function(j, n, x, d, z, y, gamma.init.set, q1, theta.init.set, epsilon, maxit, n.init, K, c, upper, lower, lower.zero, xx, xy) {

    gamma <- gamma.init.set[, j]
    out <- em_test_gamma(n, x, d, z, y, gamma, q1, theta.init.set, epsilon, maxit, n.init, K, c, upper, lower, lower.zero, xx, xy)
    out

  }, n = n, x = x, d = d, z = z, y = y, gamma.init.set = gamma.init.set, q1 = q1, theta.init.set = theta.init.set, epsilon = epsilon, maxit = maxit, n.init = n.init, K = K, c = c, upper = upper, lower = lower, lower.zero = lower.zero, xx = xx, xy = xy)

  loglik.all <- c(t(sapply(out, "[[", "loglik")))
  theta.all <- t(sapply(out, "[[", "theta"))
  diff.all <- t(sapply(out, "[[", "diff"))
  diff.vec <- apply(diff.all, 2, min)
  index <- which.max(loglik.all)

  loglik.alt <- loglik.all[index]
  theta.alt <- theta.all[index, ]

  out.null <- homogeneousMLE(x, d, y)
  loglik.null <- out.null$ll
  theta.null <- out.null$par

  em.stat <- 2*(loglik.alt - loglik.null)
  a <- list(em.stat = em.stat, theta.null = theta.null, diff.vec = diff.vec)
}

#' A function to generate a set of initial values of gamma
#' @param q2 the dimension of gamma(z)

generate_gamma <- function(c, q2, J) {

  c1 <- c[1]
  c2 <- c[2]
  c3 <- c[3]

  mat.gamma <- matrix(double(q2*J), nrow = q2)

  mat.gamma[1, ] <- runif(J, min = 0, max = c3) * (-1)^(rbinom(J, 1, 0.5))
  mat.gamma[-1, ] <- runif(J*(q2-1), min = c1, max = c2) * matrix((-1)^(rbinom(J*(q2-1), 1, 0.5)), nrow = q2 - 1)

  mat.gamma
}

#' A function to generate a set of initial values for theta
#' @param x.dim dimension of X including intercept

theta_initial <- function(x.dim, upper, lower, lower.zero, n.init) {

  alpha.init <- matrix(runif(x.dim*n.init, min = lower, max = upper), nrow = x.dim)
  beta.init <- runif(n.init, min = lower, max = upper)
  lambda.init <- runif(n.init, min = lower, max = upper)
  sigma.init <- runif(n.init, min = lower.zero, max = upper)

  theta.init <- rbind(alpha.init, beta.init, lambda.init, sigma.init)
  theta.init
}

#' A function to implement EM algorithm with gamma fixed
#' @param n sample size
#' @param theta.init the initial value of theta for EM algorithm

EM_gamma_fixed <- function(n, x, d, z, y, gamma, q1, theta.init, epsilon, maxit, xx, xy) {


  alpha <- theta.init[1:q1]
  beta <- theta.init[q1+1]
  lambda <- theta.init[q1+2]
  sigma <- theta.init[q1+3]
  diff <- 1.0
  oldloglik <- -Inf
  x_alpha <- x %*% alpha
  z_gamma <- z %*% gamma
  w.mat <- matrix(nrow=n, ncol=2)
  mean.1 <- y - x_alpha - d*beta
  mean.2 <- y - x_alpha - d*lambda

  # alpha <- theta.init[1:q1]
  # beta <- theta.init[(q1 + 1)]
  # lambda <- theta.init[(q1 + 2)]
  # sigma <- theta.init[(q1 + 3)]
  # pi <- plogis(z %*% gamma)
  # pi.alt <- 1 - pi
  # w.mat <- matrix(nrow = n, ncol = 2)

  # an initial E-step
  # mean.2 <- y - x%*%alpha - d*beta
  # mean.1 <- mean.2 - d*lambda
  # w.mat[, 1] <- pi * dnorm(mean.1, 0, sigma)
  # w.mat[, 2] <- pi.alt * dnorm(mean.2, 0, sigma)
  # l.vec <- rowSums(w.mat)
  # w.mat <- w.mat / l.vec
  # loglik <- sum(log(l.vec))
  # oldloglik <- loglik

  for(i1 in 1:maxit) {

    ##E-step
    log.prop <- plogis(z_gamma)
    w.mat[, 1] <- log.prop * dnorm(mean.1, 0, sigma)
    w.mat[, 2] <- (1-log.prop) * dnorm(mean.2, 0, sigma)
    l.vec <- rowSums(w.mat)
    w.mat <- w.mat / l.vec

    ##Exit loop or not
    loglik <- sum(log(l.vec))
    diff <- loglik - oldloglik
    oldloglik <- loglik

    #print(diff)
    if(diff < epsilon || loglik == -Inf){
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
    mat.equation.beta <- cbind(t(x_w1_d), sum(w1_d), 0)
    mat.equation.lambda <- cbind(t(x_w2_d), 0, sum(w2_d))
    mat.equation <- rbind(mat.equation.alpha, mat.equation.beta, mat.equation.lambda)

    vec.equation.alpha <- xy
    vec.equation.beta <- sum(y * w1_d)
    vec.equation.lambda <- sum(y * w2_d)
    vec.equation <- c(vec.equation.alpha, vec.equation.beta, vec.equation.lambda)

    coef <- try(solve(mat.equation, vec.equation), silent = TRUE)
    if(class(coef) == "try-error"){
      break
    }

    alpha <- coef[1:q1]
    beta <- coef[q1+1]
    lambda <- coef[q1+2]
    x_alpha <- x %*% alpha

    mean.1 <- y - x_alpha - d*beta
    mean.2 <- y - x_alpha - d*lambda

    sigma <- sqrt(mean(w1 * (mean.1)^2) + mean(w2 * (mean.2)^2))

    # # M-step
    # w <- w.mat[, 1]
    # coef_vec <- c(t(x) %*% y, sum(y * d), sum(w * y * d))
    # data_mat_w <- cbind(x, d, w * d)
    # coef_mat <- t(data_mat_w) %*% data_mat_w
    # coef_mat[x.dim + 2, x.dim + 2] <- sum(w * d)
    # const_mat <- matrix(0, nrow = x.dim + 2, ncol = x.dim + 2)
    # const_mat[x.dim + 2, x.dim + 2] <- 1
    # result <- try(solve.QP(coef_mat, coef_vec, const_mat), silent = TRUE)
    # if(class(result) == "try-error"){
    #   loglik <- -Inf
    #   break
    # }
    # #result <- tryCatch({solve.QP(coef_mat, coef_vec, const_mat)}
    # #                   , error = function(err) {coef_mat <- nearPD(coef_mat)$mat
    # #                   return(solve.QP(coef_mat, coef_vec, const_mat))})
    # alpha <- result$solution[1:x.dim]
    # beta <- result$solution[x.dim + 1]
    # lambda <- result$solution[x.dim + 2]
    # mean.2 <- y - x %*% alpha - d * beta
    # mean.1 <- mean.2 - d*lambda
    # sigma <- sqrt(mean(w * (mean.1)^2) + mean((1 - w) * (mean.2)^2))

    # # E-step
    # w.mat[, 1] <- pi * dnorm(mean.1, 0, sigma)
    # w.mat[, 2] <- pi.alt * dnorm(mean.2, 0, sigma)
    # l.vec <- rowSums(w.mat)
    # w.mat <- w.mat / l.vec
    #
    # ##Exit loop or not
    # loglik <- sum(log(l.vec))
    # diff <- loglik - oldloglik
    # oldloglik <- loglik
    #
    #
    # if(diff < epsilon){
    #   break
    # }

  }

  theta <- c(alpha, beta, lambda, sigma)
  a <- list(theta = theta, alpha = alpha, beta = beta, lambda = lambda, sigma = sigma, loglik = loglik, diff = diff)

  a

}

#' A function to implement calculate l_1 (eta^K_j) in equation (13) of Shen and He(2015)
#' @param n sample size
#' @param theta.init the initial value of theta for EM algorithm

em_test_gamma <- function(n, x, d, z, y, gamma, q1, theta.init.set, epsilon, maxit, n.init, K, c, upper, lower, lower.zero, xx, xy, n.core = detectCores() - 1) {

  #failure.theta <- 0
  c1 <- c[1]
  c2 <- c[2]
  c3 <- c[3]

  out.mle <- lapply(1:n.init, function(j, theta.init.set, n, x, d, z, y, gamma, q1, epsilon, maxit, xx, xy){

    theta.init <- theta.init.set[, j]
    EM_gamma_fixed_cpp(n, x, d, z, y, gamma, q1, theta.init, epsilon, maxit, xx, xy)

  }, theta.init.set = theta.init.set, n = n, x = x, d = d, z = z, y = y, gamma = gamma, q1 = q1, epsilon = epsilon, maxit = maxit, xx = xx, xy = xy)

  loglik.all <- c(t(sapply(out.mle, "[[", "loglik")))
  theta.all <- t(sapply(out.mle, "[[", "theta"))
  diff.all <- c(t(sapply(out.mle, "[[", "diff")))
  diff.min1 <- min(diff.all, na.rm = TRUE)
  index <- which.max(loglik.all)

  theta0 <- theta.all[index, ]
  loglik0 <- loglik.all[index]
  alpha0 <- theta0[1:q1]
  beta0 <- theta0[(q1 + 1)]
  lambda0 <- theta0[(q1 + 2)]
  sigma0 <- theta0[(q1 + 3)]
  gamma0 <- gamma

  opts <- list("algorithm"="NLOPT_LD_LBFGS", "xtol_rel"=1.0e-8)

  # an initial E-step
  w.mat <- matrix(nrow = n, ncol = 2)
  x_alpha0 <- x%*%alpha0
  mean.1 <- y - x_alpha0 - d*beta0
  mean.2 <- y - x_alpha0 - d*lambda0
  pi <- plogis(z %*% gamma0)
  w.mat[, 1] <- pi * dnorm(mean.1, 0, sigma0)
  w.mat[, 2] <- (1 - pi) * dnorm(mean.2, 0, sigma0)
  l.vec <- rowSums(w.mat)
  w.mat <- w.mat / l.vec

  for(k in 1:K){

    oldloglik0 <- loglik0
    # M-step
    w1 <- w.mat[, 1]
    w2 <- w.mat[, 2]
    w1_d <- w1 * d
    w2_d <- w2 * d
    x_w1_d <- t(x) %*% w1_d
    x_w2_d <- t(x) %*% w2_d
    mat.equation.alpha <- cbind(xx, x_w1_d, x_w2_d)
    mat.equation.beta <- cbind(t(x_w1_d), sum(w1_d), 0)
    mat.equation.lambda <- cbind(t(x_w2_d), 0, sum(w2_d))
    mat.equation <- rbind(mat.equation.alpha, mat.equation.beta, mat.equation.lambda)

    vec.equation.alpha <- xy
    vec.equation.beta <- sum(y * w1_d)
    vec.equation.lambda <- sum(y * w2_d)
    vec.equation <- c(vec.equation.alpha, vec.equation.beta, vec.equation.lambda)

    coef <- try(solve(mat.equation, vec.equation), silent = TRUE)
    if(class(coef) == "try-error"){
      break
    }

    alpha0 <- coef[1:q1]
    beta0 <- coef[q1+1]
    lambda0 <- coef[q1+2]
    x_alpha0 <- x %*% alpha0

    mean.1 <- y - x_alpha0 - d*beta0
    mean.2 <- y - x_alpha0 - d*lambda0

    sigma <- sqrt(mean(w1 * (mean.1)^2) + mean(w2 * (mean.2)^2))

    # # M-step
    # w <- w.mat[, 1]
    # coef_vec <- c(t(x) %*% y, sum(y * d), sum(w * y * d))
    # data_mat_w <- cbind(x, d, w * d)
    # coef_mat <- t(data_mat_w) %*% data_mat_w
    # coef_mat[x.dim + 2, x.dim + 2] <- sum(w * d)
    # const_mat <- matrix(0, nrow = x.dim + 2, ncol = x.dim + 2)
    # const_mat[x.dim + 2, x.dim + 2] <- 1
    # result <- try(solve.QP(coef_mat, coef_vec, const_mat), silent = TRUE)
    # if(class(result) == "try-error"){
    #   #loglik0 <- -Inf
    #   #failure.theta <- 1
    #   break
    # }
    # alpha0 <- result$solution[1:x.dim]
    # beta0 <- result$solution[x.dim + 1]
    # lambda0 <- result$solution[x.dim + 2]
    # mean.2 <- y - x %*% alpha0 - d * beta0
    # mean.1 <- mean.2 - d*lambda0
    # sigma0 <- sqrt(mean(w * (mean.1)^2) + mean((1 - w) * (mean.2)^2))

    gamma0.old <- gamma0
    out.logit <- try(nloptr::nloptr(x0 = gamma0, eval_f = llogit, eval_grad = dlogit, opts = opts, z=z, w=w1), silent = TRUE)
    if(class(out.logit) == "try-error"){
      w.mat[, 1] <- pi * dnorm(mean.1, 0, sigma0)
      w.mat[, 2] <- (1 - pi) * dnorm(mean.2, 0, sigma0)
      l.vec <- rowSums(w.mat)
      w.mat <- w.mat / l.vec
      loglik0 <- sum(log(l.vec))
      diff <- loglik0 - oldloglik0
      break
    }
    gamma0 <- out.logit$solution
    gamma0_intercept <- gamma0[1]
    gamma0_covariate <- gamma0[-1]
    max_gamma0_cov <- max(abs(gamma0_covariate))
    if(abs(gamma0_intercept) >= c3 || max_gamma0_cov >= c2 || max_gamma0_cov <= c1) {
      gamma0 <- gamma0.old
    }

    # E-step
    pi <- plogis(z %*% gamma0)
    w.mat[, 1] <- pi * dnorm(mean.1, 0, sigma0)
    w.mat[, 2] <- (1 - pi) * dnorm(mean.2, 0, sigma0)
    l.vec <- rowSums(w.mat)
    w.mat <- w.mat / l.vec
    #oldloglik0 <- loglik0
    loglik0 <- sum(log(l.vec))
    diff <- loglik0 - oldloglik0
    if(diff < 0){
      break
    }
    #print(loglik0 - oldloglik0)
  }
  diff.min2 <- diff


  theta0 <- c(alpha0, beta0, lambda0, sigma0)
  theta.init.set <- theta_initial(q1, upper, lower, lower.zero, n.init)
  theta.init.set <- cbind(theta0, theta.init.set)

  out.final.mle <- lapply(1:(n.init+1), function(j, theta.init.set, n, x, d, z, y, gamma, q1, epsilon, maxit, xx, xy){

    theta.init <- theta.init.set[, j]
    EM_gamma_fixed_cpp(n, x, d, z, y, gamma, q1, theta.init, epsilon, maxit, xx, xy)

  }, theta.init.set = theta.init.set, n = n, x = x, d = d, z = z, y = y, gamma = gamma0, q1 = q1, epsilon = epsilon, maxit = maxit, xx = xx, xy = xy)

  loglik.final.all <- c(t(sapply(out.final.mle, "[[", "loglik")))
  theta.final.all <- t(sapply(out.final.mle, "[[", "theta"))
  diff.final.all <- c(t(sapply(out.final.mle, "[[", "diff")))
  diff.min3 <- min(diff.final.all, na.rm = TRUE)
  index <- which.max(loglik.all)

  theta0 <- theta.final.all[index, ]
  loglik0 <- loglik.final.all[index]
  alpha0 <- theta0[1:q1]
  beta0 <- theta0[(q1 + 1)]
  lambda0 <- theta0[(q1 + 2)]
  sigma0 <- theta0[(q1 + 3)]

  a <- list(loglik = loglik0, theta = theta0, alpha = alpha0, beta = beta0, lambda = lambda0, sigma = sigma0, diff = c(diff.min1, diff.min2, diff.min3))

  a
}

#' A function to implement hypothesis tests with bootstrap
#' @param B the number of bootstrap replications

bootstrap_test <- function(x, d, z, y, J, c, K, B, n.init = 20, upper = 5, lower = -5, lower.zero = 0.05, epsilon = 1e-08, maxit = 500, n.core = detectCores() - 1) {

  xx <- t(x) %*% x
  xy <- t(x) %*% y

  out.em <- em_test(x, d, z, y, J, c, K, xx, xy)

  #q1 <- ncol(x)
  #q2 <- ncol(z)

  alpha.null <- out.em$theta.null$alpha
  beta.null <- out.em$theta.null$beta
  sigma.null <- out.em$theta.null$sigma
  n <- length(y)

  em.stat.boot <- mclapply(1:B, function(b, x, d, z, J, c, K, alpha.null, beta.null, sigma.null, n, xx, xy) {

    y <- x %*% alpha.null + d * beta.null + rnorm(n, 0, sigma.null)
    xy <- t(x) %*% y
    out.em.b <- em_test(x, d, z, y, J, c, K, xx, xy)
    #print(b)
    #print(out.em.b$diff.vec)
    out.em.b$em.stat


  }, x = x, d = d, z = z, J = J, c = c, K = K, alpha.null = alpha.null, beta.null = beta.null, sigma.null = sigma.null, n = n, xx = xx, xy = xy, mc.cores = n.core)

  #set.seed(1)
  #for(b in 1:238){
  #  print(b)
  #  y <- x %*% alpha.null + d * beta.null + rnorm(n, 0, sigma.null)
  #  xy <- t(x) %*% y
  #  out.em.b <- em_test(x, d, z, y, J, c, K, xx, xy)
  #}
  #em.stat.boot <- lapply(1:B, function(b, x, d, z, J, c, K, alpha.null, beta.null, sigma.null, n, xx, xy) {

  #  y <- x %*% alpha.null + d * beta.null + rnorm(n, 0, sigma.null)
  #  xy <- t(x) %*% y
  #  out.em.b <- em_test(x, d, z, y, J, c, K, xx, xy)
    #print(b)
    #print(out.em.b$diff.vec)
  #  out.em.b$em.stat


  #}, x = x, d = d, z = z, J = J, c = c, K = K, alpha.null = alpha.null, beta.null = beta.null, sigma.null = sigma.null, n = n, xx = xx, xy = xy)


  em.stat.boot <- unlist(em.stat.boot)

  em.stat <- out.em$em.stat

  p.val <- mean(em.stat < em.stat.boot)

  a <- list(em.stat = em.stat, p.val = p.val)

  a
}

