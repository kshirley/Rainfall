########################################################################################################################
# Assume that you have done the data setup, and you've fit the model using new/mcmc_baseline.R to either the 
# simulated data or the real data.
# Here we'll look at trace plots of the posterior draws

# Clear the workspace:
rm(list=ls())
setwd("~/Git/Rainfall/")

#setwd("~kathrynvasilaky/Documents/OneDrive/IRI/RainfallSimulation/Rainfall/Rainfall")
#path<-"kathrynvasilaky/Documents/OneDrive/IRI/RainfallSimulation/Rainfall/Rainfall"

#read in scripts:
source("Rcode_tobit_mcmc_functions.R")
library(MASS)
library(msm)
library(mvtnorm)
library(coda)

# load the input data and assign elements of the list to the global namespace:
load("input_data.RData")
for (i in 1:length(input.data)) assign(names(input.data)[i], input.data[[i]])

# load the simulated data and true values of parameters and assign to global namespace:
#seed <- 888
#load(paste0("sim_data_seed_", seed, ".RData"))
#for (i in 1:length(sim.data)) assign(names(sim.data)[i], sim.data[[i]])

# load the mcmc fit to this simluated data:
# This was a 2000 iteration, 3-chain MCMC to simulated data seed = 888:
#gibbs_out_20150306_G2k.RData, gibbs_out_20150616_G20k.RData
#gibbs_out_20150616_G20k.RData
#gibbs_out_20150514_G10k.RData
load(file = "gibbs_out_20150616_G20k.RData")
for (i in 1:length(gibbs.list)) assign(names(gibbs.list)[i], gibbs.list[[i]])

# set up a few parameters related to the MCMC itself:
adapt <- 2000
G <- 20000
K <- 3



### beta:
pdf(file = "plots-20k/fig_tobit_trace_beta_P23.pdf", width = 8, height = 6)
par(mfrow = c(2, 3))
for (j in 1:P){
  for (s in 1:S){
    tp(beta.gibbs[, , j, s], thin= 20, ylim = range(beta.gibbs[ , , j, ]), 
       main = paste(X.names[j], site.names[s], sep=" "), las = 1, xlab = "Iteration")
  }
}
dev.off()

# All the betas look good re: convergence

# G-R statistics:
gr.beta <- matrix(0, P, S)  # store them in a 23 x 5 matrix
for (s in 1:S) {
  x <- as.mcmc.list(list(mcmc(beta.gibbs[1, , 1:P, s], start = 10001, end = G, thin = 1), 
                         mcmc(beta.gibbs[2, , 1:P, s], start = 10001, end = G, thin = 1), 
                         mcmc(beta.gibbs[3, , 1:P, s], start = 10001, end = G, thin = 1)))
  gr.beta[, s] <- gelman.diag(x)$psrf[, 1]
}  

round(gr.beta, 2)
apply(gr.beta, 2, range)
apply(gr.beta, 1, range)


### mu
pdf(file = "plots-20k/fig_tobit_trace_mu.pdf", width = 8, height = 6)
par(mfrow = c(1, 1))
for (j in 1:P){
  tp(mu.gibbs[, , j], thin = 20, main = X.names[j], las = 1, xlab = "Iteration")
}
dev.off()

# look OK

# G-R statistics:
x <- as.mcmc.list(list(mcmc(mu.gibbs[1, , 1:P], start = 10001, end = G, thin = 1), 
                       mcmc(mu.gibbs[2, , 1:P], start = 10001, end = G, thin = 1), 
                       mcmc(mu.gibbs[3, , 1:P], start = 10001, end = G, thin = 1)))
gr.mu <- gelman.diag(x)$psrf[, 1]



### sigma
pdf(file = "plots-20k/fig_tobit_trace_sigma.pdf", width = 8, height = 6)
par(mfrow = c(1, 1))
for (j in 1:P){
  tp(sigma.gibbs[, , j], thin = 20, main = X.names[j], las = 1, xlab = "Iteration")
}
dev.off()

# G-R statistics:
x <- as.mcmc.list(list(mcmc(sigma.gibbs[1, , 1:P], start = 10001, end = G, thin = 1), 
                       mcmc(sigma.gibbs[2, , 1:P], start = 10001, end = G, thin = 1), 
                       mcmc(sigma.gibbs[3, , 1:P], start = 10001, end = G, thin = 1)))
gr.sigma <- gelman.diag(x)$psrf[, 1]




### tau:
tp(tau.gibbs, las = 1, main = "tau", ylab = "tau")
tp(tau.gibbs, las = 1, main = "tau", ylab = "tau", burn = 10000)

gelman.diag(as.mcmc.list(list(mcmc(as.matrix(tau.gibbs[1, ]), start = 10001, end = G, thin = 1), 
                              mcmc(as.matrix(tau.gibbs[2, ]), start = 10001, end = G, thin = 1), 
                              mcmc(as.matrix(tau.gibbs[3, ]), start = 10001, end = G, thin = 1))))
# 1.06                              



### lambda:

# basic trace plot:
tp(lambda.gibbs, las = 1, main = "lambda", ylab = "lambda")

# zooming in on post burnin:
tp(lambda.gibbs, las = 1, main = "lambda", ylab = "lambda", burn = 7500)

# explicitly limiting y axis to c(0, 0.2)
tp(lambda.gibbs, las = 1, main = "lambda", ylab = "lambda", ylim = c(0, 0.2))

gelman.diag(as.mcmc.list(list(mcmc(as.matrix(lambda.gibbs[1, ]), start = 10001, end = G, thin = 1), 
                              mcmc(as.matrix(lambda.gibbs[2, ]), start = 10001, end = G, thin = 1), 
                              mcmc(as.matrix(lambda.gibbs[3, ]), start = 10001, end = G, thin = 1))))
# 1.19





### mu_arc:
tp(mu.arc.gibbs, las=1, main="mu arc", ylab="mu arc")
gelman.diag(as.mcmc.list(list(mcmc(as.matrix(mu.arc.gibbs[1, ]), start = 10001, end = G, thin = 1), 
                              mcmc(as.matrix(mu.arc.gibbs[2, ]), start = 10001, end = G, thin = 1), 
                              mcmc(as.matrix(mu.arc.gibbs[3, ]), start = 10001, end = G, thin = 1))))
# 1.00

### tau_arc
tp(tau.arc.gibbs, las = 1, main = "mu arc", ylab = "tau arc")
gelman.diag(as.mcmc.list(list(mcmc(as.matrix(tau.arc.gibbs[1, ]), start = 10001, end = G, thin = 1), 
                              mcmc(as.matrix(tau.arc.gibbs[2, ]), start = 10001, end = G, thin = 1), 
                              mcmc(as.matrix(tau.arc.gibbs[3, ]), start = 10001, end = G, thin = 1))))
# 1.01


### beta.arc
pdf(file = "plots-20k/fig_tobit_trace_beta_arc.pdf", width = 8, height = 6)
par(mfrow = c(1, 1))
for (s in 1:S){
  tp(beta.arc.gibbs[, , s], thin = 20, main = site.names[s], las = 1, xlab = "Iteration")
}
dev.off()

# look OK

# G-R statistics:
x <- as.mcmc.list(list(mcmc(beta.arc.gibbs[1, , 1:S], start = 10001, end = G, thin = 1), 
                       mcmc(beta.arc.gibbs[2, , 1:S], start = 10001, end = G, thin = 1), 
                       mcmc(beta.arc.gibbs[3, , 1:S], start = 10001, end = G, thin = 1)))
gr.beta.arc <- gelman.diag(x)$psrf[, 1]

range(gr.beta.arc)
# 1.01 to 1.32




### Sigma: 
pdf(file = "plots-20k/fig_tobit_trace_sigma_covariance.pdf", width = 8, height = 8)
for (s in 1:S) {
  par(mfrow = c(J[s], J[s]))
  if (s != 6) {
    for (i in 1:J[s]) {
      for (j in 1:J[s]) {  	
        tp(Sigma.gibbs[[s]][, , i, j], las = 1, thin = 20, xlab = "Iteration")
        title(main = paste0("Site = ", s, ", Row = ", i, ", Col = ", j))
      }
    }
  } else {
    for (i in 1:J[s]) {
      for (j in 1:J[s]) {
      	par(mar = c(0, 2, 0, 0), xaxt = "n")
        tp(Sigma.gibbs[[s]][, , i, j], las = 1, thin = 20, xlab = "", main = "", ylab = "")
      }
    }  	
  }
}
dev.off()



# posterior mean of time series for spatial mean:
mu.bar <- apply(mu.gibbs[, 10001:G, ], 3, mean)

# spatial mean across locations:
z.bar <- X %*% matrix(mu.bar, ncol = 1)

# plot the spatial mean across locations by time:
plot(date.string, z.bar, type = "l")

sel <- as.numeric(substr(date.string, 1, 4)) == 2005
plot(date.string[sel], z.bar[sel], type = "l")













# posterior means of covariance matrices
lapply(Sigma.gibbs, function(x) apply(x[, 5001:10000, , ], 3:4, mean))

# posterior means of correlation matrices:
lapply(Sigma.gibbs, function(x) cov2cor(apply(x[, 5001:10000, , ], 3:4, mean)))

# posterior means of correlation matrices:
lapply(Sigma.gibbs, function(x) sqrt(diag(apply(x[, 5001:10000, , ], 3:4, mean))))

apply(Sigma.gibbs[[1]][, 5001:10000, , ], 3:4, mean)






#Ruben Gelman Stat for a few key parameters using the Ruben Gelman statistic, modified for our lists. 
#For a little more detail on where the adapted code cam from see RGelman_test.R

burn <- 15000:20000
#Lambda
list <- as.list(rep(NA, 3))
#get the data into a list of 3 lists for the chains
for (s in 1:3) list[[s]] <- lambda.gibbs[s, burn]
RG.lambda <- kvgelman.diag(list,confidence = 0.95, transform = FALSE, autoburnin = FALSE, multivariate = TRUE)

#Mu
list <- as.list(rep(NA, 3))
RG.mu <- as.list(NA)

for (s in 1:3) {
  for (j in 1:23) {
    list[[s]] <- mu.gibbs[s,,j]
    RG.mu[j] <- kvgelman.diag(list,confidence = 0.95, transform = FALSE, autoburnin = TRUE, multivariate = TRUE)
  }
}

#Tau
list <- as.list(rep(NA, 3))
for (s in 1:3) list[[s]] <- tau.gibbs[s,burn]
RG.tau <- kvgelman.diag(list,confidence = 0.95, transform = FALSE, autoburnin = FALSE, multivariate = TRUE)

#Mu Arc
list <- as.list(rep(NA, 3))
for (s in 1:3) list[[s]] <- mu.arc.gibbs[s,burn]
RG.mu.arc <- kvgelman.diag(list,confidence = 0.95, transform = FALSE, autoburnin = FALSE, multivariate = TRUE)

#Tau Arc
list <- as.list(rep(NA, 3))
for (s in 1:3) list[[s]] <- tau.arc.gibbs[s,burn]
RG.tau.arc <- kvgelman.diag(list,confidence = 0.95, transform = FALSE, autoburnin = FALSE, multivariate = TRUE)


RG.lambda
RG.mu
RG.tau
RG.mu.arc
RG.tau.arc




# posterior means of monthly el nino effects:
apply(beta.gibbs[, 5001:10000, 12:23, ], 3:4, mean)

apply(mu.gibbs[, 5001:10000, ], 3, mean)
apply(sigma.gibbs[, 5001:10000, ], 3, mean)

apply(beta.gibbs[, 5001:10000, 12:23, ], 3:4, mean)

tp(beta.gibbs[, , 8, 2])



