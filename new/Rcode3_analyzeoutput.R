########################################################################################################################

# Assume that you have done the data setup, and you've fit the model using new/mcmc_baseline.R to either the 
# simulated data or the real data.

# Here we'll look at trace plots of the posterior draws

# Clear the workspace:
rm(list=ls())
setwd("~/Git/Rainfall/")

#setwd("~/SkyDrive/IRI/RainfallSimulation/Rainfall")
#path<-"~/SkyDrive/IRI/RainfallSimulation/Rainfall/"

# read in scripts:
source("Rcode_tobit_mcmc_functions.R")
library(MASS)
library(msm)
library(mvtnorm)

# load the input data and assign them to the global namespace:
load("input_data.RData")
for (i in 1:length(input.data)) assign(names(input.data)[i], input.data[[i]])

# load the simulated data and true values of parameters and assign to global namespace:
seed <- 888
load(paste0("sim_data_seed_", seed, ".RData"))
for (i in 1:length(sim.data)) assign(names(sim.data)[i], sim.data[[i]])

# load the mcmc fit to this simluated data:
# This was a 2000 iteration, 3-chain MCMC to simulated data seed = 888:
load(file = "gibbs_out_20150306_G2k.RData")
for (i in 1:length(gibbs.list)) assign(names(gibbs.list)[i], gibbs.list[[i]])

# set up a few parameters related to the MCMC itself:
adapt <- 500
G <- 2000
K <- 3

# trace plot:
tp <- function(gibbs.input,   # a (num.chains x num.iterations) matrix
               burn = 0,
               end = dim(gibbs.input)[2], 
               nc = 3,   # number of chains
               ylim = range(gibbs.input[, (burn + 1):end]), 
               thin = 1, 
               ...) {
  if(nc == 1){
  	z <- matrix(gibbs.input,ncol=1)
  } else {
    z <- matrix(t(gibbs.input),ncol=nc)
  }
  G.local <- dim(z)[1]
  thin.seq <- seq(thin, G.local, by = thin)
  lt <- length(thin.seq)
  plot(thin.seq,z[thin.seq,1],col=1,type="l",ylim=ylim,...,xaxt="n")
  axis(1,at=seq(0,G.local,length=5))
  if (nc > 1) {
    for (i in 2:nc) lines(thin.seq, z[thin.seq, i], col = i)
  }
}

tp(lambda.gibbs, las = 1, main = "lambda", ylab = "lambda")
abline(h = sim.data$lambda, col = 4)

tp(tau.gibbs, las = 1, main = "tau", ylab = "tau")
abline(h = sim.data$tau, col = 4)

pdf(file = "fig_Sigma_check.pdf", width = 5, height = 5)
for (s in 1:S) {
  l <- sum(!upper.tri(sim.data$Sigma[[s]]))
  r <- row(sim.data$Sigma[[s]])[!upper.tri(sim.data$Sigma[[s]])]
  c <- col(sim.data$Sigma[[s]])[!upper.tri(sim.data$Sigma[[s]])]
  #print(dim(sim.data$Sigma[[s]]))
  #print(paste(r, c, l))
  for (j in 1:l) {
    tp(Sigma.gibbs[[s]][, , r[j], c[j]], las = 1, ylim = c(-4, 4))
    title(main = paste0("Site = ", s, ", Row = ", r[j], ", Col = ", c[j]))
    abline(h = sim.data$Sigma[[s]][r[j], c[j]], col = 4)
  }
}
dev.off()









