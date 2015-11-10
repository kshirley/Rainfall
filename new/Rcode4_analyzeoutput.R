###############################################################################

# Assume that you have done the data setup, and you've fit the model using 
# new/mcmc_baseline.R to either the simulated data or the real data.
# Here we'll look at trace plots of the posterior draws

# Clear the workspace:
rm(list = ls())
setwd("~/Git/Rainfall/")

#setwd("~/Documents/OneDrive/IRI/RainfallSimulation/Rainfall/Rainfall")
#path<-"~/Documents/OneDrive/IRI/RainfallSimulation/Rainfall/Rainfall"

#read in scripts:
source("Rcode_tobit_mcmc_functions.R")
library(MASS)
library(msm)
library(mvtnorm)
library(coda)

# load the input data and assign elements of the list to the global namespace:
load("input_data.RData")
for (i in 1:length(input.data)) assign(names(input.data)[i], input.data[[i]])

# load the simulated data and true values of parameters and assign to global 
# namespace:
seed <- 267
load(paste0("sim_data_seed_", seed, ".RData"))
for (i in 1:length(sim.data)) assign(names(sim.data)[i], sim.data[[i]])

# load the mcmc fit to this simluated data:
# This was a 2000 iteration, 3-chain MCMC to simulated data seed = 888:

#gibbs_out_20150904_Sim20k.RData
#gibbs_out_20150616_G20k.RData
#load(file = "gibbs_out_20150616_G20k.RData")
#load(file = "gibbs_out_2015sim267_G20k.RData")
load(file = "gibbs_out_sim_267_G20k.RData")
for (i in 1:length(gibbs.list)) assign(names(gibbs.list)[i], gibbs.list[[i]])

# set up a few parameters related to the MCMC itself:
adapt <- 2000
G <- 20000
K <- 3



### beta (seasonal effects):
pdf(file = "plots-sim267/fig_tobit_trace_beta.pdf", width = 8, height = 6)
par(mfrow = c(2, 3))
for (j in 1:P){
  for (s in 1:S){
    tp(beta.gibbs[, , j, s], thin = 100, ylim = range(beta.gibbs[ , , j, ]), 
       main = paste(X.names[j], site.names[s], sep=" "), las = 1, 
       xlab = "Iteration")
    abline(h = sim.data$beta[j, s], col = 4)
  }
}
dev.off()

### mu (mean of the betas)
pdf(file = "plots-sim267/fig_tobit_trace_mu.pdf", width = 8, height = 6)
par(mfrow = c(1, 1))
for (j in 1:P){
  tp(mu.gibbs[, , j], thin = 100, main = X.names[j], las = 1, 
     xlab = "Iteration", ylim = range(mu.gibbs))
  abline(h = sim.data$mu[j], col = 4)
}
dev.off()

### sigma (sd of betas)
pdf(file = "plots-sim267/fig_tobit_trace_sigma.pdf", width = 8, height = 6)
par(mfrow = c(1, 1))
for (j in 1:P){
  tp(sigma.gibbs[, , j], thin = 100, main = X.names[j], las = 1, 
     xlab = "Iteration", ylim = range(sigma.gibbs))
  abline(h = sim.data$sigma[j], col = 4)
}
dev.off()

### tau: (tau^2*V, scaling factor of correlation matrix to produce covariance 
###       matrix for Zt, spatial mean rainfall)
pdf(file = "plots-sim267/fig_tobit_trace_tau.pdf", width = 8, height = 6)
tp(tau.gibbs, las = 1, main = "tau", ylab = "tau", thin = 100)
abline(h = sim.data$tau, col = 4)
dev.off()

### lambda: (scaler of distance correlation matrix, V=lambda*distance)
pdf(file = "plots-sim267/fig_tobit_trace_lambda.pdf", width = 8, height = 6)
tp(lambda.gibbs, las = 1, main = "lambda", ylab = "lambda", burn = (G/2 + 1), 
   thin = 100)
abline(h = sim.data$lambda, col = 4)
dev.off()

### mu_arc:
pdf(file = "plots-sim267/fig_tobit_trace_mu_arc.pdf", width = 8, height = 6)
tp(mu.arc.gibbs, las = 1, main = "mu_arc", ylab = "mu_arc", thin = 100)
abline(h = sim.data$mu.arc, col = 4)
dev.off()

### tau_arc:
pdf(file = "plots-sim267/fig_tobit_trace_tau_arc.pdf", width = 8, height = 6)
tp(tau.arc.gibbs, las = 1, main = "tau_arc", ylab = "tau_arc", thin = 100)
abline(h = sim.data$tau.arc, col = 4)
dev.off()

### beta.arc is the arc bias on the Xarc variable for each site
pdf(file = "plots-sim267/fig_tobit_trace_beta_arc.pdf", width = 8, height = 6)
par(mfrow = c(1, 1))
for (s in 1:S){
  tp(beta.arc.gibbs[, , s], thin = 100, main = site.names[s], las = 1, 
     xlab = "Iteration", ylim = range(beta.arc.gibbs))
  abline(h = sim.data$beta.arc[s], col = 4)
}
dev.off()

### Sigma: is the big spatial variation matrix between sites)
pdf(file = "plots-sim267/fig_tobit_trace_sigma_covariance.pdf", width = 8, 
    height = 8)
for (s in 1:S) {
  par(mfrow = c(J[s], J[s]))
  if (s != 6) {
    for (i in 1:J[s]) {
      for (j in 1:J[s]) {  	
        tp(Sigma.gibbs[[s]][, , i, j], las = 1, thin = 100, xlab = "Iteration")
        abline(h=sim.data$Sigma[[s]][i, j], col = 4, lwd = 2)
        title(main = paste0("Site = ", s, ", Row = ", i, ", Col = ", j))
      }
    }
  } else {
    for (i in 1:J[s]) {
      for (j in 1:J[s]) {
      	par(mar = c(0, 2, 0, 0), xaxt = "n")
        tp(Sigma.gibbs[[s]][, , i, j], las = 1, thin = 100, xlab = "", 
           main = "", ylab = "")
      	abline(h = sim.data$Sigma[[s]][i, j], col = 4, lwd = 2)
      }
    }  	
  }
}
dev.off()




# Gelman-Rubin statistics for each batch of parameters:

# beta
gr.beta <- matrix(0, P, S)  # store them in a 23 x 6 matrix
for (s in 1:S) {
  x <- as.mcmc.list(list(mcmc(beta.gibbs[1, , 1:P, s], start = (G/2 + 1), 
                              end = G, thin = 1), 
                         mcmc(beta.gibbs[2, , 1:P, s], start = (G/2 + 1), 
                              end = G, thin = 1), 
                         mcmc(beta.gibbs[3, , 1:P, s], start = (G/2 + 1), 
                              end = G, thin = 1)))
  gr.beta[, s] <- gelman.diag(x)$psrf[, 1]
}

round(gr.beta, 2)
apply(gr.beta, 2, range)
apply(gr.beta, 1, range)

# mu
x <- as.mcmc.list(list(mcmc(mu.gibbs[1, , 1:P], start = (G/2 + 1), end = G, 
                            thin = 1), 
                       mcmc(mu.gibbs[2, , 1:P], start = (G/2 + 1), end = G, 
                            thin = 1), 
                       mcmc(mu.gibbs[3, , 1:P], start = (G/2 + 1), end = G, 
                            thin = 1)))
gr.mu <- gelman.diag(x)$psrf[, 1]

# sigma
x <- as.mcmc.list(list(mcmc(sigma.gibbs[1, , 1:P], start = (G/2 + 1), end = G, 
                            thin = 1), 
                       mcmc(sigma.gibbs[2, , 1:P], start = (G/2 + 1), end = G, 
                            thin = 1), 
                       mcmc(sigma.gibbs[3, , 1:P], start = (G/2 + 1), end = G, 
                            thin = 1)))
gr.sigma <- gelman.diag(x)$psrf[, 1]

# tau
gr.tau <- gelman.diag(as.mcmc.list(list(
  mcmc(as.matrix(tau.gibbs[1, ]), start = (G/2 + 1), end = G, thin = 1), 
  mcmc(as.matrix(tau.gibbs[2, ]), start = (G/2 + 1), end = G, thin = 1), 
  mcmc(as.matrix(tau.gibbs[3, ]), start = (G/2 + 1), end = G, thin = 1)))
)
# 1.06                              

# lambda
gr.lambda <- gelman.diag(as.mcmc.list(list(
  mcmc(as.matrix(lambda.gibbs[1, ]), start = (G/2 + 1), end = G, thin = 1), 
  mcmc(as.matrix(lambda.gibbs[2, ]), start = (G/2 + 1), end = G, thin = 1), 
  mcmc(as.matrix(lambda.gibbs[3, ]), start = (G/2 + 1), end = G, thin = 1)))
)
# 1.19

# mu_arc
gr.mu.arc <- gelman.diag(as.mcmc.list(list(
  mcmc(as.matrix(mu.arc.gibbs[1, ]), start = (G/2 + 1), end = G, thin = 1), 
  mcmc(as.matrix(mu.arc.gibbs[2, ]), start = (G/2 + 1), end = G, thin = 1), 
  mcmc(as.matrix(mu.arc.gibbs[3, ]), start = (G/2 + 1), end = G, thin = 1)))
)

# tau_arc
gr.tau.arc <- gelman.diag(as.mcmc.list(list(
  mcmc(as.matrix(tau.arc.gibbs[1, ]), start = (G/2 + 1), end = G, thin = 1), 
  mcmc(as.matrix(tau.arc.gibbs[2, ]), start = (G/2 + 1), end = G, thin = 1), 
  mcmc(as.matrix(tau.arc.gibbs[3, ]), start = (G/2 + 1), end = G, thin = 1)))
)

# beta_arc
x <- as.mcmc.list(list(
  mcmc(beta.arc.gibbs[1, , 1:S], start = (G/2 + 1), end = G, thin = 1), 
  mcmc(beta.arc.gibbs[2, , 1:S], start = (G/2 + 1), end = G, thin = 1), 
  mcmc(beta.arc.gibbs[3, , 1:S], start = (G/2 + 1), end = G, thin = 1))
)
gr.beta.arc <- gelman.diag(x)$psrf[, 1]

# first site and first list of correlations for the site over 3 chains
# Sigma.gibbs[[site]][chain, runs , correlation matrix , column/series]
gr.Sigma <- list("vector", S)
for (s in 1:S) {
  gr.Sigma.gibbs[[s]] <- matrix(0, J[s], J[s])
  for (i in 1:J[s]) {
    for (j in 1:J[s]) {    
      x <- as.mcmc.list(list(
             mcmc(Sigma.gibbs[[s]][1, ,i , j], start = (G/2 + 1), thin = 1), 
             mcmc(Sigma.gibbs[[s]][2, ,i , j], start = (G/2 + 1), thin = 1), 
             mcmc(Sigma.gibbs[[s]][3, ,i , j], start = (G/2 + 1), thin = 1))
      )
      gr.Sigma[[s]][i, j] <- gelman.diag(x)$psrf[1]
    }
  }
}







RG1 <- matrix(gr.Sigma.gibbs[1:4], nrow = 2, byrow = TRUE)
RG2 <- matrix(gr.Sigma.gibbs[5:8], nrow = 2, byrow = TRUE)
RG3 <- matrix(gr.Sigma.gibbs[9:12], nrow = 2, byrow = TRUE)
RG4 <- matrix(gr.Sigma.gibbs[13:16], nrow = 2, byrow = TRUE)
RG5 <- matrix(gr.Sigma.gibbs[17:20], nrow = 2, byrow = TRUE)
RG6 <- matrix(gr.Sigma.gibbs[21:45], nrow = 5, byrow = TRUE)



# Gelman-Rubin stats
gr.beta.table <- xtable(gr.beta)
gr.mu.table <- xtable(as.matrix(gr.mu))

mat4 <- matrix(0,5,2)
mat4[1,2] <- gr.tau[[1]][1]
mat4[1,1] <- "gr.tau"
mat4[2,2] <- gr.lambda[[1]][1]
mat4[2,1] <- "gr.lambda"
mat4[3,2] <- gr.tau.arc[[1]][1]
mat4[3,1] <- "gr.tau.arc"
mat4[4,2] <- gr.beta.arc[[1]][1]
mat4[4,1] <- "gr.beta.arc"
mat4[5,2] <- gr.sigma
mat4[5,1] <- "gr.sigma"
mat4.table <- xtable(mat4)

RG1.table <- xtable(RG1)
RG2.table <- xtable(RG2)
RG3.table <- xtable(RG3)
RG4.table <- xtable(RG4)
RG5.table <- xtable(RG5)
RG6.table <- xtable(RG6)



### Plots comparing the priors to the posteriors (for a few parameters):
post.burn <- 10001:20000
dense <- density(lambda.gibbs[, post.burn])

plot(density(rgamma(100000, shape = 1, scale = 1)), xlab = "lambda", 
     ylab = "p(lambda)", main = "lambda", xlim = c(0, 5), las = 1, 
     ylim = c(0, max(dense$y)))
abline(h = 0, col = gray(0.5))
lines(dense, col = 2)
legend("topright", inset = 0.01, col = c(1, 2), 
       legend = c("Prior", "Posterior"), lwd = 1)







mean.hat <- 10
var.hat <- 10
v0.tau <- 2*mean.hat^2/var.hat + 4  # 24
s0.tau <- mean.hat*(v0.tau - 2)/v0.tau  # 55/6

plot(density(sqrt(1/rgamma(100000, shape = v0.tau/2, scale = v0.tau*s0.tau/2))), 
     xlab = "tau", ylab = "p(tau)", main = "tau")


# posterior mean of lambda
#lambda.hat <- mean(lambda.gibbs[, (G/2 + 1):G])
#plot(seq(0, 1, by = 0.01)*100, exp(-lambda.hat*seq(0, 1, by = 0.01)), 
#     xlab = "km")






# posterior mean of time series for spatial mean:
mu.bar <- apply(mu.gibbs[, 10001:G, ], 3, mean)

# spatial mean across locations:
z.bar <- X %*% matrix(mu.bar, ncol = 1)

# plot the spatial mean across locations by time:
plot(date.string, z.bar, type = "l")
sel <- as.numeric(substr(date.string, 1, 4)) == 2005
plot(date.string[sel], z.bar[sel], type = "l")



# posterior mean of time series for spatial mean:
time.bar <- apply(mu.gibbs[, 10001:G, 1:11], 3, mean)

# spatial mean across locations:
z.time <- X[, 1:11] %*% matrix(time.bar, ncol = 1)

# plot the spatial mean across locations by time:
plot(date.string, z.time, type = "l")
sel <- as.numeric(substr(date.string, 1, 4)) == 2005
plot(date.string[sel], z.time[sel], type = "l")
abline(h = 0, lty = 2)

f.sd <- function(x) {
  sqrt(diag(apply(x[, 10001:G, , ], 3:4, mean)))
}

round(unlist(sapply(Sigma.gibbs, f.sd)), 1)














# posterior means of covariance matrices
lapply(Sigma.gibbs, function(x) apply(x[, 5001:10000, , ], 3:4, mean))

# posterior means of correlation matrices:
lapply(Sigma.gibbs, function(x) cov2cor(apply(x[, 5001:10000, , ], 3:4, mean)))

# posterior means of correlation matrices:
lapply(Sigma.gibbs, function(x) sqrt(diag(apply(x[, 5001:10000, , ], 3:4, mean))))

apply(Sigma.gibbs[[1]][, 5001:10000, , ], 3:4, mean)






# Gelman-Rubin stat for a few key parameters
# For more detail on where the adapted code came from see RGelman_test.R

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



