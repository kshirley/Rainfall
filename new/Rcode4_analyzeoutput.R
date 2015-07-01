########################################################################################################################
# Assume that you have done the data setup, and you've fit the model using new/mcmc_baseline.R to either the 
# simulated data or the real data.
# Here we'll look at trace plots of the posterior draws

# Clear the workspace:
rm(list=ls())
setwd("~/Git/Rainfall/")

setwd("~kathrynvasilaky/Documents/OneDrive/IRI/RainfallSimulation/Rainfall/Rainfall")
path<-"kathrynvasilaky/Documents/OneDrive/IRI/RainfallSimulation/Rainfall/Rainfall"

#read in scripts:
source("Rcode_tobit_mcmc_functions.R")
library(MASS)
library(msm)
library(mvtnorm)
library(coda)

# load the input data and assign them to the global namespace:
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

tp(lambda.gibbs, las = 1, main = "lambda", ylab = "lambda")
#abline(h = sim.data$lambda, col = 4)
tp(lambda.gibbs, las = 1, main = "lambda", ylab = "lambda", burn=7500)

#tau
tp(tau.gibbs, las = 1, main = "tau", ylab = "tau", burn=7500)
#abline(h = sim.data$tau, col = 4)

#Mu arc
tp(mu.arc.gibbs, las=1, main="mu arc",ylab="mu arc")

#Tau arc
tp(tau.arc.gibbs, las=1, main="mu arc",ylab="tau arc")

#Sigma 
pdf(file = "fig_Sigma_check.pdf", width = 5, height = 5)
for (s in 1:S) {
  l <- sum(!upper.tri(Sigma.gibbs[[s]][1, 1, , ]))
  r <- row(Sigma.gibbs[[s]][1, 1, , ])[!upper.tri(Sigma.gibbs[[s]][1, 1, , ])]
  c <- col(Sigma.gibbs[[s]][1, 1, , ])[!upper.tri(Sigma.gibbs[[s]][1, 1, , ])]
  #print(dim(sim.data$Sigma[[s]]))
  #print(paste(r, c, l))
  for (j in 1:l) {
    tp(Sigma.gibbs[[s]][, , r[j], c[j]], las = 1)
    title(main = paste0("Site = ", s, ", Row = ", r[j], ", Col = ", c[j]))
    #abline(h = sim.data$Sigma[[s]][r[j], c[j]], col = 4)
  }
}
dev.off()

# posterior means of covariance matrices
lapply(Sigma.gibbs, function(x) apply(x[, 5001:10000, , ], 3:4, mean))

# posterior means of correlation matrices?
lapply(Sigma.gibbs, function(x) cov2cor(apply(x[, 5001:10000, , ], 3:4, mean)))

# posterior means of correlation matrices?
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



