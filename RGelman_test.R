#http://www.sumsar.net/blog/2013/06/three-ways-to-run-bayesian-models-in-r/
#http://www.johnmyleswhite.com/notebook/2010/08/29/mcmc-diagnostics-in-r-with-the-coda-package/
#http://dr-k-lo.blogspot.com/2013/03/problems-loading-jags-into-r.html


setwd("~kathrynvasilaky/Documents/OneDrive/new_afterlostmac")

library(rjags)
library(coda)
library(plyr)

#Create data of 100 points
N <- 10

#Create two lists
x <- as.list(rep(NA,2))

#Create 3 chains
#change 1 to 2 for more than one variable
x[[1]] <- matrix(c(rnorm(20,0,1)), N, 1)
x[[2]] <- matrix(c(rnorm(5,0,1)), N, 1)
x[[3]] <- matrix(c(rnorm(40,0,1)), N, 1)

#Change the starting points
x[[1]][1] <- 30
x[[1]][2] <- 50

x[[2]][1] <- -30
x[[2]][2] <- -50

x[[3]][1] <- .05
x[[3]][2] <- -.05


#Write my own function 
kvgelman.diag <- function (x, confidence = 0.95, transform = FALSE, autoburnin = TRUE, 
          multivariate = TRUE) 
{
  #x <- as.mcmc.list(x)
  if (length(x) < 2) 
    stop("You need at least two chains")
  #if (autoburnin && start(x) < end(x)/2) 
   # x <- window(x, start = end(x)/2 + 1)
  
  Niter <- length(x[[1]])
  #Niter <- niter(x)
  Nchain <- length(x)
  #Nchain <- nchain(x)
  Nvar <- length(x[[1]][1])
  #Nvar <- nvar(x)
  
  xnames <- varnames(x)
  #if (transform) 
   # x <- gelman.transform(x)
  #put x into a list---this is where my function can start with a regular list
  x <- lapply(x, as.matrix)
  #sapply generates the variance and covariance between the means of each chain
  #array places them into 3 matrices
  S2 <- array(sapply(x, var, simplify = TRUE), dim = c(Nvar, 
                                                       Nvar, Nchain))
  W <- apply(S2, c(1, 2), mean)
  xbar <- matrix(sapply(x, apply, 2, mean, simplify = TRUE), 
                 nrow = Nvar, ncol = Nchain)
  B <- Niter * var(t(xbar))
  if (Nvar > 1 && multivariate) {
    if (is.R()) {
      CW <- chol(W)
      emax <- eigen(backsolve(CW, t(backsolve(CW, B, transpose = TRUE)), 
                              transpose = TRUE), symmetric = TRUE, only.values = TRUE)$values[1]
    }
    else {
      emax <- eigen(qr.solve(W, B), symmetric = FALSE, 
                    only.values = TRUE)$values
    }
    mpsrf <- sqrt((1 - 1/Niter) + (1 + 1/Nvar) * emax/Niter)
  }
  else mpsrf <- NULL
  w <- diag(W)
  b <- diag(B)
  s2 <- matrix(apply(S2, 3, diag), nrow = Nvar, ncol = Nchain)
  muhat <- apply(xbar, 1, mean)
  var.w <- apply(s2, 1, var)/Nchain
  var.b <- (2 * b^2)/(Nchain - 1)
  cov.wb <- (Niter/Nchain) * diag(var(t(s2), t(xbar^2)) - 2 * 
                                    muhat * var(t(s2), t(xbar)))
  V <- (Niter - 1) * w/Niter + (1 + 1/Nchain) * b/Niter
  var.V <- ((Niter - 1)^2 * var.w + (1 + 1/Nchain)^2 * var.b + 
              2 * (Niter - 1) * (1 + 1/Nchain) * cov.wb)/Niter^2
  df.V <- (2 * V^2)/var.V
  df.adj <- (df.V + 3)/(df.V + 1)
  B.df <- Nchain - 1
  W.df <- (2 * w^2)/var.w
  R2.fixed <- (Niter - 1)/Niter
  R2.random <- (1 + 1/Nchain) * (1/Niter) * (b/w)
  R2.estimate <- R2.fixed + R2.random
  R2.upper <- R2.fixed + qf((1 + confidence)/2, B.df, W.df) * 
    R2.random
  psrf <- cbind(sqrt(df.adj * R2.estimate), sqrt(df.adj * R2.upper))
  dimnames(psrf) <- list(xnames, c("Point est.", "Upper C.I."))
  out <- list(psrf = psrf, mpsrf = mpsrf)
  class(out) <- "gelman.diag"
  out
}


kvgelman.diag(x,confidence = 0.95, transform = FALSE, autoburnin = TRUE, 
            multivariate = TRUE)


#########################################################################

#k <- 3
#x <- as.list(rep(NA,k))

#for (i in 1:k) {
#  x[[i]] <- rnorm(10, 0, 1)
#}


N <- 100
x <- 1:N
epsilon <- rnorm(N, 0, 1)
y <- x + epsilon


# The model specification
model_string <- "model{
  for(i in 1:length(y)) {
    y[i] ~ dnorm(mu, tau)
  }
  mu ~ dnorm(0, 0.0001)
  sigma ~ dlnorm(0, 0.0625)
  tau <- 1 / pow(sigma, 2)
}"

# Running the model
model <- jags.model(textConnection(model_string), data = list(y = y), n.chains = 3, n.adapt= 10000)
update(model, 100); # Burnin for 10000 samples
mcmc_samples <- coda.samples(model, variable.names=c("mu", "sigma"), n.iter=200)
x <- mcmc_samples






gelman.diag <- function (x, confidence = 0.95, transform = FALSE, autoburnin = TRUE, 
                           multivariate = TRUE) 
{
  x <- as.mcmc.list(x)
  if (nchain(x) < 2) 
    stop("You need at least two chains")
  if (autoburnin && start(x) < end(x)/2) 
    x <- window(x, start = end(x)/2 + 1)
  
  #Niter <- length(xtest[[1]][,1])
  Niter <- niter(x)
  #Nchain <- length(xtest)
  Nchain <- nchain(x)
  
  #Nvar <- length(xtest[[1]][1,])
  Nvar <- nvar(x)
  
  xnames <- varnames(x)
  if (transform) 
    x <- gelman.transform(x)
  #put x into a list---this is where my function can start with a regular list
  x <- lapply(x, as.matrix)
  #sapply generates the variance and covariance between the means of each chain
  #array places them into 3 matrices
  S2 <- array(sapply(x, var, simplify = TRUE), dim = c(Nvar, 
                                                       Nvar, Nchain))
  W <- apply(S2, c(1, 2), mean)
  xbar <- matrix(sapply(x, apply, 2, mean, simplify = TRUE), 
                 nrow = Nvar, ncol = Nchain)
  B <- Niter * var(t(xbar))
  if (Nvar > 1 && multivariate) {
    if (is.R()) {
      CW <- chol(W)
      emax <- eigen(backsolve(CW, t(backsolve(CW, B, transpose = TRUE)), 
                              transpose = TRUE), symmetric = TRUE, only.values = TRUE)$values[1]
    }
    else {
      emax <- eigen(qr.solve(W, B), symmetric = FALSE, 
                    only.values = TRUE)$values
    }
    mpsrf <- sqrt((1 - 1/Niter) + (1 + 1/Nvar) * emax/Niter)
  }
  else mpsrf <- NULL
  w <- diag(W)
  b <- diag(B)
  s2 <- matrix(apply(S2, 3, diag), nrow = Nvar, ncol = Nchain)
  muhat <- apply(xbar, 1, mean)
  var.w <- apply(s2, 1, var)/Nchain
  var.b <- (2 * b^2)/(Nchain - 1)
  cov.wb <- (Niter/Nchain) * diag(var(t(s2), t(xbar^2)) - 2 * 
                                    muhat * var(t(s2), t(xbar)))
  V <- (Niter - 1) * w/Niter + (1 + 1/Nchain) * b/Niter
  var.V <- ((Niter - 1)^2 * var.w + (1 + 1/Nchain)^2 * var.b + 
              2 * (Niter - 1) * (1 + 1/Nchain) * cov.wb)/Niter^2
  df.V <- (2 * V^2)/var.V
  df.adj <- (df.V + 3)/(df.V + 1)
  B.df <- Nchain - 1
  W.df <- (2 * w^2)/var.w
  R2.fixed <- (Niter - 1)/Niter
  R2.random <- (1 + 1/Nchain) * (1/Niter) * (b/w)
  R2.estimate <- R2.fixed + R2.random
  R2.upper <- R2.fixed + qf((1 + confidence)/2, B.df, W.df) * 
    R2.random
  psrf <- cbind(sqrt(df.adj * R2.estimate), sqrt(df.adj * R2.upper))
  dimnames(psrf) <- list(xnames, c("Point est.", "Upper C.I."))
  out <- list(psrf = psrf, mpsrf = mpsrf)
  class(out) <- "gelman.diag"
  out
}



gelman.diag(x,confidence = 0.95, transform = FALSE, autoburnin = TRUE, 
            multivariate = TRUE)




###


list <- as.list(rep(NA, 3))
for (s in 1:3) list[[s]] <- mu.gibbs[s,,22]

kvgelman.diag(list,confidence = 0.95, transform = FALSE, autoburnin = TRUE, multivariate = TRUE)


