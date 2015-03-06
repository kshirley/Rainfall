########################################################################################################################

# Start here after running 'Rcode1_data_setup.R'
rm(list=ls())
setwd("~/Git/Rainfall/")

# load the input data and assign them to the global namespace:
load("input_data.RData")
for (i in 1:length(input.data)) assign(names(input.data)[i], input.data[[i]])

# read in scripts:
source("Rcode_tobit_mcmc_functions.R")
library(MASS)

# Simulate data:

# Set true parameter values:
set.seed(836)
mu <- rnorm(P, 0, 1)
sigma <- runif(P, 0.3, 0.8)
alpha <- 5 # degrees of freedom for student-t
lambda <- 3 # spatial correlation decay
tau <- 4 # scalar for spatial covariance
mu.arc <- 1 # mean ARC effect
beta.arc <- seq(0, 2, length = S)
tau.arc <- sd(beta.arc) # variability of ARC effects

# Simulate regression coefficients for underlying weather process:
set.seed(328)
beta <- matrix(NA, P, S)
for (p in 1:P){
  beta[p, ] <- rnorm(S, mu[p], sigma[p])
}

# Set up covariance matrices for series within each site:
Sigma <- as.list(rep(NA,S))
Sigma[[1]] <- matrix(c(1, 0.8, 0.8, 1), 2, 2)
Sigma[[2]] <- matrix(c(1, 0.4, 0.4, 1), 2, 2)
Sigma[[3]] <- matrix(c(3, 1, 1, 3), 2, 2)
Sigma[[4]] <- matrix(c(1, -0.6, -0.6, 1), 2, 2)
Sigma[[5]] <- diag(J[5])
Sigma[[6]] <- diag(J[6])

# create exponential covariance matrix:
V <- R.cov(lambda, d.mat)

# Simulate Z given X_t, beta, and V (lambda and tau)
Z <- mvrnorm(n = N, mu = rep(0, S), Sigma = tau^2*V)
xbeta <- X %*% beta
Z <- Z + xbeta

# Draw random multivariate student-t variables:
W <- as.list(rep(NA, S))
gamma.mat <- matrix(NA, N, S)
for (s in 1:S){
  #print(s)
  gamma <- rgamma(n = N, shape = alpha/2, scale = 2/alpha)
  W[[s]] <- mvrnorm(n = N, mu = rep(0, J[s]), Sigma = Sigma[[s]])/sqrt(gamma)
  W[[s]] <- t(W[[s]] + matrix(Z[, s], N, J[s]) + matrix(rep(X.arc[[s]]*beta.arc[s], each = N), ncol = J[s]))
  gamma.mat[, s] <- gamma
}

# Truncate at zero, and delete values to match missing data properties of the real data:
Y.complete <- as.list(rep(NA, S))  # with no missing data, but truncated
Y.sim <- as.list(rep(NA, S))  # with some data missing
for (s in 1:S){
  Y.complete[[s]] <- W[[s]]
  Y.complete[[s]][Y.complete[[s]] < 0] <- 0
  Y.sim[[s]] <- Y.complete[[s]]
  Y.sim[[s]][na.mat[[s]] == 1] <- NA
}
# Great. Now Y.sim is like our real data

# save the true parameter values in a list:
sim.data <- list(mu = mu,
                 sigma = sigma,
                 alpha = alpha,
                 lambda = lambda,
                 tau = tau,
                 mu.arc = mu.arc,
                 tau.arc = tau.arc,
                 beta.arc = beta.arc,
                 beta = beta,
                 Sigma = Sigma,
                 Z = Z,
                 W = W,
                 gamma.mat = gamma.mat, 
                 Y.sim = Y.sim)

save(sim.data, file = "sim_data_seed_836.RData")


