library(dplyr)
library(mvtnorm)
library(polynom)
library(ggplot2)
library(MASS)

women <- read.table("https://raw.githubusercontent.com/STIMALiU/BayesLearnCourse/master/Labs/WomenWork.dat", header=TRUE)
model <- glm(Work ~ 0 + ., data = women, family = binomial) # fit logistic regression
#model
target <- as.vector(women[, 1])
covariates <- as.matrix(women[, 2:9])
nPara <- ncol(covariates)
tau <- 10
mu0 <- as.vector(rep(0,nPara))
sigma0 <- tau^2*diag(nPara)

logprior <- function(betas,sigma0) {
  mu <- matrix(rep(0,8),ncol = 1)
  prior <- dmvnorm(betas, mean=mu, sigma=sigma0, log=TRUE)
  return(prior)
}

loglik <- function(betas, target, covariates) {
  LP <- covariates%*%betas
  proba <- LP*target - log(1 + exp(LP))
  logliklihood <- sum(proba)
  if (abs(logliklihood) == Inf) {
    logliklihood = -20000
  }
  return(logliklihood)
}

logpost<-function(betas, target, covariates, mu0, sigma0) {
  betas <- as.vector(betas)
  logposter <- logprior(betas, sigma0) + loglik(betas, target, covariates)
  return(logposter)
}

initVal <- rep(0,8)
OptimResults <- optim(initVal, logpost,
                      target = target, covariates = covariates, mu0 = mu0, sigma0 = sigma0,
                      method="BFGS", control=list(fnscale=-1), hessian=TRUE)

OptimResults$par # the optimal values of each parameter

## Simulate from predictive distribution

data <- matrix(c(1, 10, 8, 10, (10/10)^2, 40, 1, 1),ncol = 8)
sigma <- solve(OptimResults$hessian)
mu <- OptimResults$par
prdictions<-function() {
  draws<-c()
  for(i in 1:1000) {
    betas<-rmvnorm(1,mu,sigma)
    LP<- data*betas
    prob<-exp(LP)/(1+exp(LP))
    draws[i]<-rbinom(1,1,prob)
  }
  draws
}
table(prdictions())
