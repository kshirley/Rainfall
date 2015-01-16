# Comes from R_code_tobit_real3.R, which comes from R_code_tobit_real3.R, which was split into
# real and simulated data from R code tobit.R

# Reorganizing this file so that it goes from
# 1. Data organization
#   a. Create date lists
#   b. Create X Arc matrix of weather cosines and sines
#   c. Set starting values for gibbs sampling
#   d. Set priors
#   e. Set up storage for sampling
# 2. Gibbs sampling (3 hours for 3 chains)
# 3. Save burn in data 

######################
# 1.a Read in the data
######################
rm(list=ls())

# set working directory:
#setwd("/Users/kathrynvasilaky/SkyDrive/IRI/RainfallSimulation/Rainfall")
#setwd("/Users/kshirley/Stats/Rainfall")
setwd("~/Git/Rainfall/")

# enter in here wherever you want to store the scripts and data files
#path <- getwd()#"~/SkyDrive/IRI/RainfallSimulation/Rainfall"  
#path <- paste(path,'/', sep='')
path <- "~/Git/Rainfall/"

# a lot of gibbs sampling functions:
source(paste(path, "R_code_multisite_covariance_scripts.R", sep=""))

# load libraries:
library(MASS)
library(msm)
#library(LaplacesDemon)
library(mvtnorm)

# define student t density function:
dst <- function(x, nu, mu, sigma) gamma((nu+1)/2)/(gamma(nu/2)*sqrt(nu*pi*sigma^2))*(1+1/nu*((x-mu)/sigma)^2)^(-1*(nu+1)/2)

# Load in 15 time series for Ethiopia:
load(paste(path, "ethiopia_full_data.RData", sep=""))

# quick check on number of years of observed data per time series:
count.na <- function(x) sum(is.na(x))
round((dim(data)[2] - apply(data, 1, count.na))/365.25, 1)

# Re order the sites:
data <- data[c(6,1,7,2,8,3,9,4,10,5,11,12,13,14,15),]

# Set up number of time series:
T <- dim(data)[2]

# membership of series withing sites:
site.mat <- cbind(c(1,1,2,2,3,3,4,4,5,5,6,6,6,6,6),c(1,2,1,2,1,2,1,2,1,2,1,2,3,4,5))

# ARC indicator variable for each series:
arc <- c(1,0,1,0,1,0,1,0,1,0,1,0,0,0,0)

# number of locations
S <- 6

# number of series for each location
L <- c(2,2,2,2,2,5)
L.sum <- sum(L)
month.days <- c(31,28,31,30,31,30,31,31,30,31,30,31)
month <- c(rep(rep(1:12,month.days),49),rep(1:12,month.days)[1:209])

#set the data to 1992 onward, so get rid of 31 years:
data <- data[,(365*31+1):T] # 1992 onward, all 15 sites:

# number of days
T <- dim(data)[2]

# x is an indicator of a wet day:
x <- matrix(as.numeric(data > 0), L.sum, T)

# rename x to y, and take transpose
y <- t(x)

# create the datestring beginning at Jan 1 1992:
date.string <- as.Date(1:(T+5)-1,origin="1992-01-01")
leap.year.index <- 60 + c(0,365*4+1,365*8+2,365*12+3,365*16+4)
date.string <- date.string[-leap.year.index]
month.names <- unique(months(date.string))
year.names<-unique(format(date.string,"%Y"))

month.vec <- match(months(date.string),month.names)
year.vec <- match(year.names,year.names)

# read in lat-long data:
ll <- read.table(paste(path, "lat_long_details.csv", sep=""), sep=",", header=TRUE)
ll <- ll[c(3,5,4,1,2,6), ]
d.mat <- matrix(NA,S,S)
for (i in 1:S){
  for (j in 1:S){
    d.mat[i,j] <- sqrt((ll[i,"Lat"]-ll[j,"Lat"])^2+(ll[i,"Long"]-ll[j,"Long"])^2)*108
  }
}
rownames(d.mat) <- ll[,1]
colnames(d.mat) <- ll[,1]
d.mat <- d.mat/100

# get location names:
site.names <- unlist(strsplit(rownames(data),"_")[seq(1,by=2,length=6)])[seq(1,by=2,length=6)]

### Include El Nino stuff:
# Read in Nino 3.4 index:
# KV - this index may come from here: http://www.cgd.ucar.edu/cas/catalog/climind/TNI_N34/
# data from 1871 to present, V1 is the year
nino <- read.table(paste(path, "nino34.long.data", sep=""), colClasses=rep("numeric", 13))

# vectorize data into one column
nino <- matrix(t(as.matrix(nino[,-1])),ncol=1)

# Make all dates, 12 months from 1871 to 2010
nino.dates <- paste(rep(1871:2010,each=12),unique(months(date.string)),sep=" ")

# cut dates and data off
nino <- nino[1:which(nino.dates=="2010 April")]; nino.dates <- nino.dates[1:which(nino.dates=="2010 April")]
# nino only goes through April 2010.


# subselect just the months that align with the rainfall data:
# what index should the data start at?
w <- which(nino.dates == "1991 October") # going 3 months before first rainfall month
nino <- nino[w:length(nino.dates)]
nino.dates <- nino.dates[w:1672]

# We want to impute 3 more months of nino values to stretch to the end of the daily rainfall time series:
#nino.impute <- as.numeric(predict(ar(nino), n.ahead = 3)$pred)
#nino <- c(nino, nino.impute)
#nino.dates <- c(nino.dates, paste("2010",month.names[5:7]))
# OK, now nino is length 226, starting in October 1991 (3 months before rainfall data starts) and ending in July 2010

# center nino at its mean
nino <- nino - mean(nino)
month.vec <- match(months(date.string), month.names)
month.mat <- matrix(0, T, 12)
# repeat month data 18 times, add in 6 more months
# this: c(rep(month.days,18),month.days[1:6],28), constructs 18 repeats of the months, 
# plus the first 6 months plus one more element of 28
# so nino[1:223] are the values for each of 223 months, which is the latter expression, 
# so each value of nino will be repeated the number of days in the 223 months

#Align Nino start date with 3 month lag to data's start date
month.mat[cbind(1:T, month.vec)] <- rep(nino, c(rep(month.days, 18), month.days[1:6], 28))


# Meta-parameters:
T <- 6779
S <- 6
J <- c(2, 2, 2, 2, 2, 5)
P <- 5

# index the NA elements of the data set
# KV - he is creating an array of ones and zeros somehow related to 
# the NAs in the original data set, data
# KV - so na.mat is a list of 6 lists ... the first 5 lists will be 
# 6779x2 and the 6th will be 6779x5 based on groupings of the site data
na.mat <- as.list(rep(NA,S))  # 1 = missing
for (s in 1:S){
  na.mat[[s]] <- matrix(0, J[s], T)
  na.mat[[s]][is.na(data[site.mat[, 1] == s, ])] <- 1
}

#########################
# 1.b Enter ARC indicator
########################
# KV-ARC data is http://journals.ametsoc.org/doi/abs/10.1175/JCLI-D-12-00206.1
# KV-"(ARC) project aims to create an independent climate data record of sea surface 
#temperatures (SSTs) covering recent decades that can be used for climate change analysis"
X.arc <- as.list(rep(NA, S))
for (s in 1:S){
  X.arc[[s]] <- matrix(0, J[s], 1)
  X.arc[[s]][1, 1] <- 1
}

# Real data:
# KV - A list of lists ... the same dimension as na.mat
# KV-R is the actual data from each series and each location 
R <- as.list(rep(NA, S))
for (s in 1:S) R[[s]] <- data[1:J[s] + c(0, cumsum(J))[s], ]

# Create periodic predictor variables:
# KV - for each column (10 all together) he is creating a sine 
# and a cosine
# In one year, the wave goes up and down once, plot(sin.mat[,1])
M <- 10
m.seq <- c(1:10)
sin.mat <- matrix(NA, T, M)
cos.mat <- matrix(NA, T, M)
for (m in 1:M){
  sin.mat[, m] <- sin(2*pi*(1:T)*m.seq[m]/365)
  cos.mat[, m] <- cos(2*pi*(1:T)*m.seq[m]/365)
}

# Create X_t, which is the same across sites:
# KV - might be building a regression matrix
X <- cbind(rep(1,T), (1:T)/T - 0.5, ((1:T)/T - 0.5)^2, sin.mat[, 1:4], cos.mat[, 1:4], month.mat)

# KV-P is the number of trends in the time trend variable
P <- dim(X)[2]
X.names <- c("Intercept", "Time (linear)", "Time (quadratic)", "sin-1", "sin-2", "sin-3", "sin-4",
             "cos-1", "cos-2", "cos-3", "cos-4", paste(month.names, "Nino", sep="-"))

# compute upper bounds for each observation:
up <- as.list(rep(NA,S))
for (s in 1:S){
  up[[s]] <- matrix(NA, J[s], T)
  up[[s]][na.mat[[s]] == 1] <- Inf
  up[[s]][R[[s]] == 0] <- 0
}

# Compute true/false of whether a given observation is missing or dry:
# KV - w.draw contains true if missing data or 0 value and has same shape
# as R and na.mat
# n.draw simply counts the total number of NA's or zeros
w.draw <- as.list(rep(NA, S))
n.draw <- as.list(rep(NA, S))
for (s in 1:S){
  w.draw[[s]] <- matrix(NA, J[s], T)
  n.draw[[s]] <- numeric(J[s])
  for (j in 1:J[s]){
    w.draw[[s]][j, ] <- is.na(R[[s]][j, ]) | R[[s]][j,] == 0
    n.draw[[s]][j] <- sum(R[[s]][j, ] == 0, na.rm = TRUE) + sum(is.na(R[[s]][j, ]))
  }
}


###########################
# 1.c Set starting values # 
###########################

# Set starting values for 3 MCMC chains
K <- 3

# There is where Sigma.start is loaded in!!!
# Or, load up the starting points from a previous Gibbs sampler:
load(file=paste(path, "start_list_v2.RData", sep=""))
for (i in 1:length(start.list)) assign(names(start.list)[i],start.list[[i]])

# adjust for new predictors, including quadratic term, new sine and cosine terms, and nino terms:
mu.old <- mu.start
mu.start <- matrix(NA, K, P)
mu.start[,c(1, 2, 4, 5, 8, 9)] <- mu.old
mu.start[,c(3, 6, 7, 10:23)] <- 0

sigma.old <- sigma.start
sigma.start <- matrix(NA, K, P)
sigma.start[, c(1, 2, 4, 5, 8, 9)] <- sigma.old
sigma.start[, c(3, 6, 7, 10:23)] <- 0.5

# change alpha.start:
alpha.start <- rep(10, K)
tau.start <- 8:10
beta.start <- array(NA, dim=c(K, P, S))
for (s in 1:S) beta.start[,,s] <- mu.start

# compute the mean of Z.start
xb.start <- array(NA, dim=c(K, T, S))
for (k in 1:K) xb.start[k,,] <- X %*% beta.start[k,,]

# simulate Z.start:
Z.start <- array(NA, dim=c(K, T, S))
for (k in 1:K) Z.start[k, , ] <- mvrnorm(T, rep(0, S), tau.start[k]^2*R.cov(lambda.start[k], d.mat)) + xb.start[k, , ]

# simulate gamma.start
gamma.start <- array(NA, dim=c(K, T, S))
for (k in 1:K) gamma.start[k, , ] <- matrix(rgamma(T*S, shape = alpha.start[k]/2, scale = 2/alpha.start[k]), T, S)


##################
### 1.d Set Priors
##################
v.0 <- J
Lambda.0 <- as.list(rep(NA,S))
Lambda.0.inv <- as.list(rep(NA,S))
for (s in 1:S){
  Lambda.0[[s]] <- diag(J[s])
  Lambda.0.inv[[s]] <- solve(Lambda.0[[s]])
}

# Metropolis adaptation multipliers:
A1 <- 1.1; B1 <- 1.1^(-44/56)
xtx <- t(X)%*%X
ones <- matrix(0, S*P, P)
for (s in 1:S) ones[(s-1)*P + 1:P, ] <- diag(P)
xtx.ones <- ones %*% xtx %*% t(ones)
vec <- numeric(S^2*P^2)
for (s in 1:S) vec[(s-1)*(S*P^2)+(1:(S*P^2))] <- rep(rep((s-1)*S+1:S, each=P), P)
gamma.temp <- matrix(0, T, S)
Sigma.null <- as.list(rep(NA, S))
W.null <- as.list(rep(NA, S))
for (s in 1:S) W.null[[s]] <- matrix(0, J[s], T)


###################################
#1.e Set up storage for gibbs samples:
##################################
#chains
K <- 3
#runs
G <- 20
adapt <- 500
mu.gibbs <- array(NA,dim=c(K,G,P))
sigma.gibbs <- array(NA,dim=c(K,G,P))
alpha.gibbs <- array(NA,dim=c(K,G))
lambda.gibbs <- array(NA,dim=c(K,G))
tau.gibbs <- array(NA,dim=c(K,G))
beta.gibbs <- array(NA,dim=c(K,G,P,S))
mu.arc.gibbs <- array(NA,dim=c(K,G))
tau.arc.gibbs <- array(NA,dim=c(K,G))
beta.arc.gibbs <- array(NA,dim=c(K,G,S))
Sigma.gibbs <- as.list(rep(NA,S))
for (s in 1:S) Sigma.gibbs[[s]] <- array(NA,dim=c(K,G,J[s],J[s]))


# Save samples of W and Z, if desired:
n.samp <- 50
w.samp <- sort(sample(1:T,n.samp))
Z.gibbs <- array(NA,dim=c(K,G,n.samp,S))
W.gibbs <- as.list(rep(NA,S))
for (s in 1:S) W.gibbs[[s]] <- array(NA,dim=c(K,G,J[s],n.samp))


#################################
# 2.Start the MCMC, Gibbs Sampling
################################

dyn.load(paste(path,"zdraw.so",sep=""))
dyn.load(paste(path,"draw_gamma.so",sep=""))
is.loaded("zdraw")
is.loaded("draw_gamma")

#ng <- G

set.seed(7394)
t1 <- Sys.time()
for (k in 1:K){
  print(k)
  
  # Set jump sd for RW-MH parameters:
  jump.alpha <- 0.1; jump.lambda <- 0.05; jump.sigma <- rep(0.5,S); jump.ta <- 0.35
  
  # Fill in starting values:
  mu <- mu.start[k,]; sigma <- sigma.start[k,]; alpha <- alpha.start[k]
  lambda <- lambda.start[k]; tau <- tau.start[k]; beta <- beta.start[k,,]
  beta.arc <- beta.arc.start[k,]; mu.arc <- mu.arc.start[k]; tau.arc <- tau.arc.start[k];
  Z <- Z.start[k,,]; gamma.mat <- gamma.start[k,,]
  Sigma <- as.list(rep(NA,S))
  for (s in 1:S) Sigma[[s]] <- Sigma.start[[s]][k,,]
  
  # Store parameter values:
  mu.gibbs[k,1,] <- mu; sigma.gibbs[k,1,] <- sigma; alpha.gibbs[k,1] <- alpha
  lambda.gibbs[k,1] <- lambda; tau.gibbs[k,1] <- tau; beta.gibbs[k,1,,] <- beta
  beta.arc.gibbs[k,1,] <- beta.arc; mu.arc.gibbs[k,1] <- mu.arc; tau.arc.gibbs[k,1] <- tau.arc
  for (s in 1:S) Sigma.gibbs[[s]][k,1,,] <- Sigma[[s]]  
  Z.gibbs[k,1,,] <- Z[w.samp,]
  
  # Last, set W.start inside the loop for each chain:
  W <- draw.W.start(R,Z,beta.arc,Sigma,gamma.mat,X.arc,up,w.draw,W.null)
  for (s in 1:S) W.gibbs[[s]][k,1,,] <- W[[s]][,w.samp]
  
  # Sample W | Sigma, mu (Z, X.arc, beta.arc), gamma:
  #W <- W.start
  #for (s in 1:S) W.gibbs[[s]][k,1,,] <- W[[s]][,w.samp]
  
  # Loop through iterations:
  for (g in 2:G){
    if (g%%100==0) print(g)
    
    # Sample gamma | W, mu (Z, X.arc, beta.arc), Sigma, alpha:    
    gamma.mat <- draw.gamma(S,Z,X.arc,beta.arc,T,J,Sigma,W,alpha,gamma.temp=matrix(0,T,S))
    
    # Draw Z | W, gamma, beta.arc, Sigma, beta, lambda, tau:
    Z <- draw.Z(S,Sigma,T,J,X.arc,beta.arc,W,lambda,d.mat,tau,X,beta,gamma.mat)
    Z.gibbs[k,g,,] <- Z[w.samp,]
    
    # Draw beta.arc | W, gamma, Z, Sigma, mu.arc, tau.arc:
    beta.arc <- draw.beta.arc(S,Sigma,gamma.mat,X.arc,Z,T,J,W,tau.arc,mu.arc)
    beta.arc.gibbs[k,g,] <- beta.arc
    
    # draw mu.arc | beta.arc, tau.arc
    mu.arc <- rnorm(1,(sum(beta.arc)/tau.arc^2)/(1 + S/tau.arc^2),sqrt(1/(1 + S/tau.arc^2)))
    mu.arc.gibbs[k,g] <- mu.arc
    
    # draw tau.alpha | alpha, mu.alpha
    temp <- ta.draw(beta.arc,mu.arc,tau.arc,jump.ta,g,adapt)
    tau.arc <- temp$tau.alpha; jump.ta <- temp$jump.ta
    tau.arc.gibbs[k,g] <- tau.arc
    
    # Sample Sigma | W, mu (Z, X.arc, beta.arc), gamma:
    Sigma <- draw.Sigma(S,Z,X.arc,beta.arc,T,J,W,gamma.mat,v.0,Lambda.0.inv,Sigma.null)
    for (s in 1:S) Sigma.gibbs[[s]][k,g,,] <- Sigma[[s]]
    
    # Sample alpha | gamma.mat, a, b where prior(alpha) ~ gamma(shape=a,scale=b):
    #temp <- alpha.draw(gamma.mat,a=15,b=1,jump.alpha)
    #temp <- alpha.draw.unif(gamma.mat,a=2,b=30,jump.alpha)
    #alpha <- temp$alpha
    #jump.alpha <- temp$jump.alpha
    alpha.gibbs[k,g] <- alpha
    
    # draw beta | w, lambda:
    temp <- draw.beta.tobit(tau,lambda,vec,xtx.ones,Z,X,sigma,mu)
    beta <- temp$beta
    beta.gibbs[k,g,,] <- beta
    
    # draw tau | Z, beta, lambda, a, b where prior(tau^2) ~ inverse-gamma(shape=a,rate=b)
    tau <- tau.draw.tobit(Z,beta,lambda,a=1,b=1)
    tau.gibbs[k,g] <- tau
    
    # draw lambda | beta, tau, Z:
    temp <- lambda.draw.tobit(X,beta,tau,lambda,d.mat,Z,S,jump.lambda)
    lambda <- temp$lambda
    jump.lambda <- temp$jump.lambda
    lambda.gibbs[k,g] <- lambda
    
    # draw mu, sigma | beta:
    temp <- mu.sigma.draw(beta,mu,sigma,S,P,jump.sigma,g,adapt)
    jump.sigma <- temp$jump.sigma
    mu <- temp$mu; mu.gibbs[k,g,] <- mu
    sigma <- temp$sigma; sigma.gibbs[k,g,] <- sigma
    
    # Sample W | Sigma, mu (Z, X.arc, beta.arc), gamma:
    W <- draw.W(R,Z,beta.arc,Sigma,gamma.mat,X.arc,W,up,w.draw)
    for (s in 1:S) W.gibbs[[s]][k,g,,] <- W[[s]][,w.samp]
  }
}
t2 <- Sys.time()
t2-t1





gibbs.list <- list(mu.gibbs=mu.gibbs,sigma.gibbs=sigma.gibbs,alpha.gibbs=alpha.gibbs,lambda.gibbs=lambda.gibbs,tau.gibbs=tau.gibbs,beta.gibbs=beta.gibbs,
                   mu.arc.gibbs=mu.arc.gibbs,tau.arc.gibbs=tau.arc.gibbs,beta.arc.gibbs=beta.arc.gibbs,Sigma.gibbs=Sigma.gibbs,W.gibbs=W.gibbs)
save(gibbs.list,file=paste(path,"gibbs_out_04272014_G5000.RData",sep=""))

#load(file=paste(path,"gibbs_out_NA11172014_G5000.RData",sep=""))
load(file="gibbs_out_04272014_G5000.RData")
for (i in 1:length(gibbs.list)) assign(names(gibbs.list)[i],gibbs.list[[i]])

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

