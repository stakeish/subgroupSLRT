rm(list=ls())

library(subgroupSLRT)
library(parallel)
library(abind)
library(Matrix)
library(quadprog)
library(nloptr)
library(lbfgs)
library(speff2trial)
load("coefficients.RData")

data("ACTG175")
intercept <- coefficients[1]
coef <- coefficients[2]

data <- ACTG175[ACTG175$arms == 1 | ACTG175$arms == 3, ]
x <- cbind(rep(1, nrow(data)), data$age, data$wtkg, data$karnof, data$cd40, data$cd80, data$hemo, data$homo, data$drugs, data$race, data$gender, data$str2, data$symptom)
z <- x
dim_z <- ncol(z)
for(i1 in 2:(dim_z-1)){
  z1 <- z[, i1]
  for(i2 in (i1+1):dim_z){
    z2 <- z[, i2]
    z <- cbind(z, z1*z2)
    #print(c(i1, i2))
  }
}
d <- as.numeric(data$arms == 1)
y <- data$cd420
z <- apply(z, 2, scale)
z[, 1] <- 1

n <- nrow(data)
gamma.dim <- ncol(z)
p <- intercept + coef * n^(7/8) * sqrt(log(gamma.dim))

critical_value <- qnorm(0.95)^2

set.seed(10)

out <- SLRT_unparallel(x, d, z, y, p, ninits = 10, ninits.equal = 10)
SLRT <- out$SLRT
pval <- 1-pnorm(sqrt(pmax(0,SLRT)))

save(pval, file = "pval.RData")
