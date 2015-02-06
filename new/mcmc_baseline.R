########################################################################################################################

# Run 'Rcode1_data_setup.R'
# Now you have saved a local copy of "input_data.RData" which contains the real data and some other objects

# Next, run 'Rcode2_simulate_data.R' to simulate some data:
# Now you have another .RData file containing true parameter values and some simulated data

# Clear the workspace:
rm(list=ls())
setwd("~/Git/Rainfall/")

# read in scripts:
source("Rcode_tobit_mcmc_functions.R")
library(MASS)
library(msm)
library(mvtnorm)

# load the input data and assign them to the global namespace:
load("input_data.RData")
for (i in 1:length(input.data)) assign(names(input.data)[i], input.data[[i]])

# load the simulated data and true values of parameters and assign to global namespace:
load("sim_data_seed_836.RData")
for (i in 1:length(sim.data)) assign(names(sim.data)[i], sim.data[[i]])

# Replace the real data, Y, with the simulated version, Y.sim:
Y <- Y.sim

# Now let's write the MCMC code:
# Eventually this block might be its own file, or function called "Rcode3"

# compute upper bounds for each observation:
up <- as.list(rep(NA, S))
for (s in 1:S){
  up[[s]] <- matrix(NA, J[s], N)
  up[[s]][na.mat[[s]] == 1] <- Inf
  up[[s]][Y[[s]] == 0] <- 0
}

# Compute true/false of whether a given observation is missing or dry:
# KV - w.draw contains true if missing data or 0 value and has same shape
# as R and na.mat
# n.draw simply counts the total number of NA's or zeros
w.draw <- as.list(rep(NA, S))
n.draw <- as.list(rep(NA, S))
for (s in 1:S){
  w.draw[[s]] <- matrix(NA, J[s], N)
  n.draw[[s]] <- numeric(J[s])
  for (j in 1:J[s]){
    w.draw[[s]][j, ] <- is.na(Y[[s]][j, ]) | Y[[s]][j, ] == 0
    n.draw[[s]][j] <- sum(Y[[s]][j, ] == 0, na.rm = TRUE) + sum(is.na(Y[[s]][j, ]))
  }
}




###########################
# 1.c Set starting values # 
###########################

# Set starting values for 3 MCMC chains
K <- 3

# check names of start value parameters:
#load("start_list_v2.RData")

# random starting points:
set.seed(634)
mu.start <- matrix(rnorm(K*P), K, P)
sigma.start <- matrix(runif(K*P, 0.2, 1), K, P)
alpha.start <- rep(5, K)
lambda.start <- c(0.5, 2.5, 4.5)
tau.start <- c(2, 4, 6)
beta.start <- array(NA, dim=c(K, P, S))
for (s in 1:S) beta.start[, , s] <- mu.start + rnorm(K*P, 0, 0.5*sigma.start)
mu.arc.start <- numeric(K)
tau.arc.start <- c(0.5, 1, 1.5)
beta.arc.start <- matrix(0, K, S)
Sigma.start <- vector(mode = "list", length = S)
for (s in 1:S) {
  Sigma.start[[s]] <- array(0, dim=c(K, J[s], J[s]))
  for (k in 1:K) {
    Sigma.start[[s]][k, , ] <- diag(J[s])
  }
}


# compute the mean of Z.start
xb.start <- array(NA, dim = c(K, N, S))
for (k in 1:K) xb.start[k, , ] <- X %*% beta.start[k, , ]

# simulate Z.start:
Z.start <- array(NA, dim=c(K, N, S))
for (k in 1:K) Z.start[k, , ] <- mvrnorm(N, rep(0, S), tau.start[k]^2*R.cov(lambda.start[k], d.mat)) + xb.start[k, , ]

# simulate gamma.start
gamma.start <- array(NA, dim=c(K, N, S))
for (k in 1:K) gamma.start[k, , ] <- matrix(rgamma(N*S, shape = alpha.start[k]/2, scale = 2/alpha.start[k]), N, S)


##################
### 1.d Set Priors
##################
v.0 <- J
Lambda.0 <- as.list(rep(NA, S))
Lambda.0.inv <- as.list(rep(NA, S))
for (s in 1:S){
  Lambda.0[[s]] <- diag(J[s])
  Lambda.0.inv[[s]] <- solve(Lambda.0[[s]])
}

# Metropolis adaptation multipliers:
A1 <- 1.1; B1 <- 1.1^(-44/56)
xtx <- t(X) %*% X
ones <- matrix(0, S*P, P)
for (s in 1:S) ones[(s-1)*P + 1:P, ] <- diag(P)
xtx.ones <- ones %*% xtx %*% t(ones)
vec <- numeric(S^2*P^2)
for (s in 1:S) vec[(s - 1)*(S*P^2) + (1:(S*P^2))] <- rep(rep((s - 1)*S + 1:S, each = P), P)
gamma.temp <- matrix(0, N, S)
Sigma.null <- as.list(rep(NA, S))
W.null <- as.list(rep(NA, S))
for (s in 1:S) W.null[[s]] <- matrix(0, J[s], N)


###################################
#1.e Set up storage for gibbs samples:
##################################

K <- 3  # number of chains
G <- 100  # number of iterations/samples
adapt <- 500
mu.gibbs <- array(NA, dim=c(K, G, P))
sigma.gibbs <- array(NA, dim=c(K, G, P))
alpha.gibbs <- array(NA, dim=c(K, G))
lambda.gibbs <- array(NA, dim=c(K, G))
tau.gibbs <- array(NA, dim=c(K, G))
beta.gibbs <- array(NA, dim=c(K, G, P, S))
mu.arc.gibbs <- array(NA, dim=c(K, G))
tau.arc.gibbs <- array(NA, dim=c(K, G))
beta.arc.gibbs <- array(NA, dim=c(K, G, S))
Sigma.gibbs <- as.list(rep(NA, S))
for (s in 1:S) Sigma.gibbs[[s]] <- array(NA, dim=c(K, G, J[s], J[s]))


# Save samples of W and Z, if desired:
n.samp <- 50
w.samp <- sort(sample(1:N, n.samp))
Z.gibbs <- array(NA, dim=c(K, G, n.samp, S))
W.gibbs <- as.list(rep(NA, S))
for (s in 1:S) W.gibbs[[s]] <- array(NA, dim=c(K, G, J[s], n.samp))


#################################
# 2.Start the MCMC, Gibbs Sampling
################################

dyn.load("zdraw.so")
dyn.load("draw_gamma.so")
is.loaded("zdraw")
is.loaded("draw_gamma")

set.seed(7394)
t1 <- Sys.time()
for (k in 1:K){
  print(k)
  
  # Set jump sd for RW-MH parameters:
  jump.alpha <- 0.1
  jump.lambda <- 0.05
  jump.sigma <- rep(0.5, S)
  jump.ta <- 0.35
  
  # Fill in starting values:
  mu <- mu.start[k, ]
  sigma <- sigma.start[k, ]
  alpha <- alpha.start[k]
  lambda <- lambda.start[k]
  tau <- tau.start[k]
  beta <- beta.start[k, , ]
  beta.arc <- beta.arc.start[k, ]
  mu.arc <- mu.arc.start[k]
  tau.arc <- tau.arc.start[k]
  Z <- Z.start[k, , ]
  gamma.mat <- gamma.start[k, , ]
  Sigma <- as.list(rep(NA, S))
  for (s in 1:S) Sigma[[s]] <- Sigma.start[[s]][k, , ]
  
  # Store parameter values:
  mu.gibbs[k, 1, ] <- mu
  sigma.gibbs[k, 1, ] <- sigma
  alpha.gibbs[k, 1] <- alpha
  lambda.gibbs[k, 1] <- lambda
  tau.gibbs[k, 1] <- tau
  beta.gibbs[k, 1, ,] <- beta
  beta.arc.gibbs[k, 1, ] <- beta.arc
  mu.arc.gibbs[k, 1] <- mu.arc
  tau.arc.gibbs[k, 1] <- tau.arc
  for (s in 1:S) Sigma.gibbs[[s]][k, 1, ,] <- Sigma[[s]]  
  Z.gibbs[k, 1, , ] <- Z[w.samp, ]
  
  # Last, set W.start inside the loop for each chain:
  W <- draw.W.start(R = Y, Z, beta.arc, Sigma, gamma.mat, X.arc, up, w.draw, W.null)
  for (s in 1:S) W.gibbs[[s]][k, 1, , ] <- W[[s]][, w.samp]
  
  # Sample W | Sigma, mu (Z, X.arc, beta.arc), gamma:
  #W <- W.start
  #for (s in 1:S) W.gibbs[[s]][k,1,,] <- W[[s]][,w.samp]
  
  # Loop through iterations:
  for (g in 2:G){
    if (g %% 10 == 0) print(g)
    
    # Sample gamma | W, mu (Z, X.arc, beta.arc), Sigma, alpha:    
    gamma.mat <- draw.gamma(S, Z, X.arc, beta.arc, N, J, Sigma, W, alpha, gamma.temp=matrix(0, N, S))
    
    # Draw Z | W, gamma, beta.arc, Sigma, beta, lambda, tau:
    Z <- draw.Z(S, Sigma, N, J, X.arc, beta.arc, W, lambda, d.mat, tau, X, beta, gamma.mat)
    Z.gibbs[k, g, , ] <- Z[w.samp, ]
    
    # Draw beta.arc | W, gamma, Z, Sigma, mu.arc, tau.arc:
    beta.arc <- draw.beta.arc(S, Sigma, gamma.mat, X.arc, Z, N, J, W, tau.arc, mu.arc)
    beta.arc.gibbs[k, g, ] <- beta.arc
    
    # draw mu.arc | beta.arc, tau.arc
    mu.arc <- rnorm(1, (sum(beta.arc)/tau.arc^2)/(1 + S/tau.arc^2), sqrt(1/(1 + S/tau.arc^2)))
    mu.arc.gibbs[k, g] <- mu.arc
    
    # draw tau.alpha | alpha, mu.alpha
    temp <- ta.draw(beta.arc, mu.arc, tau.arc, jump.ta, g, adapt)
    tau.arc <- temp$tau.alpha
    jump.ta <- temp$jump.ta
    tau.arc.gibbs[k, g] <- tau.arc
    
    # Sample Sigma | W, mu (Z, X.arc, beta.arc), gamma:
    Sigma <- draw.Sigma(S, Z, X.arc, beta.arc, N, J, W, gamma.mat, v.0, Lambda.0.inv, Sigma.null)
    for (s in 1:S) Sigma.gibbs[[s]][k, g, , ] <- Sigma[[s]]
    
    # Sample alpha | gamma.mat, a, b where prior(alpha) ~ gamma(shape=a,scale=b):
    #temp <- alpha.draw(gamma.mat,a=15,b=1,jump.alpha)
    #temp <- alpha.draw.unif(gamma.mat,a=2,b=30,jump.alpha)
    #alpha <- temp$alpha
    #jump.alpha <- temp$jump.alpha
    alpha.gibbs[k, g] <- alpha
    
    # draw beta | w, lambda:
    temp <- draw.beta.tobit(tau, lambda, vec, xtx.ones, Z, X, sigma, mu)
    beta <- temp$beta
    beta.gibbs[k, g, , ] <- beta
    
    # draw tau | Z, beta, lambda, a, b where prior(tau^2) ~ inverse-gamma(shape=a,rate=b)
    tau <- tau.draw.tobit(Z, beta, lambda, a = 1, b = 1)
    tau.gibbs[k, g] <- tau
    
    # draw lambda | beta, tau, Z:
    temp <- lambda.draw.tobit(X, beta, tau, lambda, d.mat, Z, S, jump.lambda)
    lambda <- temp$lambda
    jump.lambda <- temp$jump.lambda
    lambda.gibbs[k, g] <- lambda
    
    # draw mu, sigma | beta:
    temp <- mu.sigma.draw(beta, mu, sigma, S, P, jump.sigma, g, adapt)
    jump.sigma <- temp$jump.sigma
    mu <- temp$mu
    mu.gibbs[k, g, ] <- mu
    sigma <- temp$sigma
    sigma.gibbs[k, g, ] <- sigma
    
    # Sample W | Sigma, mu (Z, X.arc, beta.arc), gamma:
    W <- draw.W(R = Y, Z, beta.arc, Sigma, gamma.mat, X.arc, W, up, w.draw)
    for (s in 1:S) W.gibbs[[s]][k, g, , ] <- W[[s]][, w.samp]
  }
}
t2 <- Sys.time()
t2 - t1


# Save the Gibbs samples in a list:
gibbs.list <- list(mu.gibbs = mu.gibbs,
                   sigma.gibbs = sigma.gibbs, 
                   alpha.gibbs = alpha.gibbs, 
                   lambda.gibbs = lambda.gibbs, 
                   tau.gibbs = tau.gibbs, 
                   beta.gibbs = beta.gibbs, 
                   mu.arc.gibbs = mu.arc.gibbs, 
                   tau.arc.gibbs = tau.arc.gibbs, 
                   beta.arc.gibbs = beta.arc.gibbs, 
                   Sigma.gibbs = Sigma.gibbs, 
                   W.gibbs = W.gibbs)
save(gibbs.list, file = "gibbs_out_20150204_G10k.RData")








#load(file=paste(path,"gibbs_out_NA11172014_G5000.RData",sep=""))
load(file = "gibbs_out_04272014_G5000.RData")
for (i in 1:length(gibbs.list)) assign(names(gibbs.list)[i], gibbs.list[[i]])




##########################
#3 New Gibbs starting points
##########################
# New starting points AFTER BURN IN:
post.burn <- 3001:G
names(gibbs.list)
mu.start <- apply(mu.gibbs[,post.burn,],c(1,3),mean)
sigma.start <- apply(sigma.gibbs[,post.burn,],c(1,3),mean)
alpha.start <- apply(alpha.gibbs[,post.burn],1,mean)
lambda.start <- apply(lambda.gibbs[,post.burn],1,mean)
tau.start <- apply(tau.gibbs[,post.burn],1,mean)
beta.start <- apply(beta.gibbs[,post.burn,,],c(1,3,4),mean)
mu.arc.start <- apply(mu.arc.gibbs[,post.burn],1,mean)
tau.arc.start <- apply(tau.arc.gibbs[,post.burn],1,mean)
beta.arc.start <- apply(beta.arc.gibbs[,post.burn,],c(1,3),mean)
Sigma.start <- as.list(rep(NA,S))
for (s in 1:S) Sigma.start[[s]] <- apply(Sigma.gibbs[[s]][,post.burn,,],c(1,3,4),mean)

start.list <- list(mu.start=mu.start,sigma.start=sigma.start,alpha.start=alpha.start,lambda.start=lambda.start,
                   tau.start=tau.start,beta.start=beta.start,mu.arc.start=mu.arc.start,tau.arc.start=tau.arc.start,
                   beta.arc.start=beta.arc.start,Sigma.start=Sigma.start)
save(start.list,file=paste(path,"start_list_v5.RData",sep=""))
# v3 is the G=5000 output from 3/12 that includes P=23 predictors and alpha=10

