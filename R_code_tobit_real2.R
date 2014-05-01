
#Reorganizing this file so that it goes from
'''
1. Data organization
  a. Create date lists
  b. Create X Arc matrix of weather cosines and sines
  c. Set starting values for gibbs sampling
  d. Set priors
  e. Set up storage for sampling
2. Gibbs sampling (3 hours for 3 chains)
3. Save burn in data 
4. Simulate new data
5. Plot historical data
6. Plot parameters from Sampling
7. Plot simulated (posterior) rainfall data
'''

######################
# 1.a Read in the data
######################
rm(list=ls())
#source("~/Stats/Misc/myRfunctions.R")
#data.path <- "~/Stats/IndexInsurance/AdiHa supporting data/"
setwd("/Users/kathrynvasilaky/SkyDrive/IRI/RainfallSimulation/R/Rainfall")
#setwd("/Users/katyav/SkyDrive/IRI/RainfallSimulation/R")
path <- getwd()#"~/SkyDrive/IRI/RainfallSimulation/R/Rainfall"  # enter in here wherever you want to store the scripts and data files
path <- paste(path,'/', sep='')
source(paste(path,"R code multisite covariance scripts.R",sep=""))  # read in some scripts I wrote
library(MASS)
library(msm)
library(LaplacesDemon)
library(mvtnorm)

# Load in old data set with Adi Ha short data sets and a couple others:
load(paste(path,"ethiopia_full_data.RData",sep=""))  # read in the data, saved as an R object
#KV-15 rows or data from 15 different sites, reorders them in 5's?
data <- data[c(6,1,7,2,8,3,9,4,10,5,11,12,13,14,15),]
T <- dim(data)[2]
site.mat <- cbind(c(1,1,2,2,3,3,4,4,5,5,6,6,6,6,6),c(1,2,1,2,1,2,1,2,1,2,1,2,3,4,5))
arc <- c(1,0,1,0,1,0,1,0,1,0,1,0,0,0,0)

#number of locations
S <- 6
#number of series for each location
L <- c(2,2,2,2,2,5)
L.sum <- sum(L)
month.days <- c(31,28,31,30,31,30,31,31,30,31,30,31)
month <- c(rep(rep(1:12,month.days),49),rep(1:12,month.days)[1:209])

# 1992 onward:
data <- data[,(365*31+1):T] # 1992 onward, all 15 sites:
#number of days
T <- dim(data)[2]
x <- matrix(as.numeric(data>0),L.sum,T)
y <- t(x)
date.string <- as.Date(1:(T+5)-1,origin="1992-01-01")
leap.year.index <- 60 + c(0,365*4+1,365*8+2,365*12+3,365*16+4)
date.string <- date.string[-leap.year.index]
month.names <- unique(months(date.string))

# read in lat-long data:
ll <- read.table(paste(path,"lat_long_details.csv",sep=""),sep=",",header=TRUE)
ll <- ll[c(3,5,4,1,2,6),]
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
nino <- read.table(paste(path,"nino34.long.data",sep=""),colClasses=rep("numeric",13))
#vectorize data into one column
nino <- matrix(t(as.matrix(nino[,-1])),ncol=1)
#Make all dates
nino.dates <- paste(rep(1871:2010,each=12),unique(months(date.string)),sep=" ")
#cut dates and data off
nino <- nino[1:1672]; nino.dates <- nino.dates[1:1672]  # nino only goes through April 2010.

# subselect just the months that align with the rainfall data:
#what index should the data start at?
w <- which(nino.dates=="1991 October") # going 3 months before first rainfall month
nino <- nino[w:length(nino.dates)]
nino.dates <- nino.dates[w:1672]

# We want to impute 3 more months of nino values to stretch to the end of the daily rainfall time series:
nino.impute <- as.numeric(predict(ar(nino),n.ahead=3)$pred)
nino <- c(nino,nino.impute)
nino.dates <- c(nino.dates,paste("2010",month.names[5:7]))
# OK, now nino is length 226, starting in October 1991 (3 months before rainfall data starts) and ending in July 2010

# center nino at its mean
nino <- nino - mean(nino)
month.vec <- match(months(date.string),month.names)
month.mat <- matrix(0,T,12)
#repeat month data 18 times, add in 6 more months
#this: c(rep(month.days,18),month.days[1:6],28), constructs 18 repeats of the months, plus the first 6 months plus one more element of 28
#so nino[1:223] are the values for each of 223 months, which is the latter expression, so each value of nino will be repeated the number of days in the 223 months

month.mat[cbind(1:T,month.vec)] <- rep(nino[1:223],c(rep(month.days,18),month.days[1:6],28))

########  Do this if you're simulating data ################################################################


########  End here if you're simulating data ################################################################
# Meta-parameters:
T <- 6779
S <- 6
J <- c(2,2,2,2,2,5)
P <- 5

# index the NA elements of the data set
# KV - he is creating an array of ones and zeros somehow related to 
# the NAs in the original data set, data
# KV - so na.mat is a list of 6 lists ... the first 5 lists will be 
# 6779x2 and the 6th will be 6779x5 based on groupings of the site data
na.mat <- as.list(rep(NA,S))  # 1 = missing
for (s in 1:S){
  na.mat[[s]] <- matrix(0,J[s],T)
  na.mat[[s]][is.na(data[site.mat[,1]==s,])] <- 1
}

#########################
# 1.b Enter ARC indicator
########################
# KV-ARC data is http://journals.ametsoc.org/doi/abs/10.1175/JCLI-D-12-00206.1
# KV-"(ARC) project aims to create an independent climate data record of sea surface 
#temperatures (SSTs) covering recent decades that can be used for climate change analysis"
X.arc <- as.list(rep(NA,S))
for (s in 1:S){
  X.arc[[s]] <- matrix(0,J[s],1)
  X.arc[[s]][1,1] <- 1
}

# Real data:
#KV - A list of lists ... the same dimension as na.mat
#KV-R is the actual data from each series and each location 
R <- as.list(rep(NA,S))
for (s in 1:S) R[[s]] <- data[1:J[s]+c(0,cumsum(J))[s],]

# Create periodic predictor variables:
# KV - for each column (10 all together) he is creating a sine 
# and a cosine
# In one year, the wave goes up and down once, plot(sin.mat[,1])
M <- 10
m.seq <- c(1:10)
sin.mat <- matrix(NA,T,M)
cos.mat <- matrix(NA,T,M)
for (m in 1:M){
  sin.mat[,m] <- sin(2*pi*(1:T)*m.seq[m]/365)
  cos.mat[,m] <- cos(2*pi*(1:T)*m.seq[m]/365)
}

# Create X_t, which is the same across sites:
# KV - might be building a regression matrix
X <- cbind(rep(1,T),(1:T)/T-0.5,((1:T)/T-0.5)^2,sin.mat[,1:4],cos.mat[,1:4],month.mat)
#KV-P is the number of trends in the time trend variable
P <- dim(X)[2]
X.names <- c("Intercept","Time (linear)","Time (quadratic)","sin-1","sin-2","sin-3","sin-4","cos-1","cos-2","cos-3","cos-4",paste(month.names,"Nino",sep="-"))

# compute upper bounds for each observation:
up <- as.list(rep(NA,S))
for (s in 1:S){
  up[[s]] <- matrix(NA,J[s],T)
  up[[s]][na.mat[[s]]==1] <- Inf
  up[[s]][R[[s]]==0] <- 0
}

# Compute true/false of whether a given observation is missing or dry:
# KV - w.draw contains true if missing data or 0 value and has same shape
# as R and na.mat
# n.draw simply counts the total number of NA's or zeros

w.draw <- as.list(rep(NA,S))
n.draw <- as.list(rep(NA,S))
for (s in 1:S){
  w.draw[[s]] <- matrix(NA,J[s],T)
  n.draw[[s]] <- numeric(J[s])
  for (j in 1:J[s]){
    w.draw[[s]][j,] <- is.na(R[[s]][j,]) | R[[s]][j,]==0
    n.draw[[s]][j] <- sum(R[[s]][j,]==0,na.rm=TRUE) + sum(is.na(R[[s]][j,]))
  }
}


######## skip this if using real data ###########
######## resume here if using real data ###########

########################
#1.c Set starting values 
########################

# Set starting values for 3 MCMC chains
K <- 3

#There is where Sigma.start is loaded in!!!
# Or, load up the starting points from a previous Gibbs sampler:
load(file=paste(path,"start_list_v2.RData",sep=""))
for (i in 1:length(start.list)) assign(names(start.list)[i],start.list[[i]])

# adjust for new predictors, including quadratic term, new sine and cosine terms, and nino terms:
mu.old <- mu.start
mu.start <- matrix(NA,K,P)
mu.start[,c(1,2,4,5,8,9)] <- mu.old
mu.start[,c(3,6,7,10:23)] <- 0

sigma.old <- sigma.start
sigma.start <- matrix(NA,K,P)
sigma.start[,c(1,2,4,5,8,9)] <- sigma.old
sigma.start[,c(3,6,7,10:23)] <- 0.5

# change alpha.start:
alpha.start <- rep(10,K)
tau.start <- 8:10
beta.start <- array(NA,dim=c(K,P,S))
for (s in 1:S) beta.start[,,s] <- mu.start

# compute the mean of Z.start
xb.start <- array(NA,dim=c(K,T,S))
for (k in 1:K) xb.start[k,,] <- X %*% beta.start[k,,]

# simulate Z.start:
Z.start <- array(NA,dim=c(K,T,S))
for (k in 1:K) Z.start[k,,] <- mvrnorm(T,rep(0,S),tau.start[k]^2*R.cov(lambda.start[k],d.mat)) + xb.start[k,,]

# simulate gamma.start
gamma.start <- array(NA,dim=c(K,T,S))
for (k in 1:K) gamma.start[k,,] <- matrix(rgamma(T*S,shape=alpha.start[k]/2,scale=2/alpha.start[k]),T,S)


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
ones <- matrix(0,S*P,P)
for (s in 1:S) ones[(s-1)*P+1:P,] <- diag(P)
xtx.ones <- ones %*% xtx %*% t(ones)
vec <- numeric(S^2*P^2)
for (s in 1:S) vec[(s-1)*(S*P^2)+(1:(S*P^2))] <- rep(rep((s-1)*S+1:S,each=P),P)
gamma.temp <- matrix(0,T,S)
Sigma.null <- as.list(rep(NA,S))
W.null <- as.list(rep(NA,S))
for (s in 1:S) W.null[[s]] <- matrix(0,J[s],T)


###################################
#1.e Set up storage for gibbs samples:
##################################
#chains
K <- 3
#runs
G <- 5000
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
# 2.Start the MCMC, Gibss Sampling
################################
source(paste(path,"R code multisite covariance scripts.R",sep=""))
dyn.load(paste(path,"zdraw.so",sep=""))
dyn.load(paste(path,"draw_gamma.so",sep=""))
is.loaded("zdraw")
is.loaded("draw_gamma")

#ng <- G

set.seed(Sys.time())
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

load(file=paste(path,"gibbs_out_04272014_G5000.RData",sep=""))
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

###############################
#4.Simulate W.new from posterior
###############################
# function to simulate rainfall given parameters from gibbs output:
sim.W <- function(alpha,beta,lambda,tau,beta.arc,Sigma,X,X.arc,na.preserve=TRUE,na.mat){
  Z.sim <- mvrnorm(T,rep(0,S),tau^2*R.cov(lambda,d.mat)) + X %*% beta
  gamma.mat <- matrix(rgamma(T*S,shape=alpha/2,scale=2/alpha),T,S)  
  W.new <- as.list(rep(NA,S))
  for (s in 1:S){
    mu.W <- t(Z.sim[,s] + t(matrix(X.arc[[s]]*beta.arc[s],J[s],T)))
    temp <- t(mvrnorm(n=T,mu=rep(0,J[s]),Sigma=Sigma[[s]])/sqrt(gamma.mat[,s]))
    temp <- temp + mu.W
    temp[temp<0] <- 0
    W.new[[s]] <- temp
    if (na.preserve) W.new[[s]][na.mat[[s]]==1] <- NA
  }
  #return(unlist(W.new)[unlist(na.mat)==0])
  return(W.new)
}

#Generate W.new data, from posterior distributions
#KV-K is number of chains, G is number of runs per chain, J.sum is the total number of sereis across locations
#KV-right now it generates 5000 runs
p.wet.gibbs <- array(NA,dim=c(K,G,J.sum,12)) # monthly percentage of wet days:
mean.gibbs <- array(NA,dim=c(K,G,J.sum,12)) # monthly percentage of wet days:
sd.gibbs <- array(NA,dim=c(K,G,J.sum,12)) # monthly percentage of wet days:
max.gibbs <- array(NA,dim=c(K,G,J.sum,12)) # monthly percentage of wet days:

zero_count_data.gibbs <- array(NA,dim = c(K,100,J.sum,6779))
for (k in 1:K){
  for (g in 1:G){
    if (g%%20==0) print(g)
    # Collect Sigma into a list:
    Sigma.sim <- as.list(rep(NA,S))
    for (s in 1:S) Sigma.sim[[s]] <- Sigma.gibbs[[s]][k,g,,]
    # Simulate rainfall:
    W.new <- sim.W(alpha=alpha.gibbs[k,g],beta=beta.gibbs[k,g,,],lambda=lambda.gibbs[k,g],tau=tau.gibbs[k,g],beta.arc=beta.arc.gibbs[k,g,],
                   Sigma=Sigma.sim,X=X,X.arc=X.arc,na.preserve=TRUE,na.mat=na.mat)
    
    if (g>4900) {
      for (s in 1:S){
        for (j in 1:J[s]){
          zero_count_data.gibbs[k,g-4900,c(0,cumsum(J))[s]+j,] <- W.new[[s]][j,]
        }
      }
    }
    
    for (s in 1:S){
      for (j in 1:J[s]){
        
        for (m in 1:12){
          #print(c(0,cumsum(J))[s]+j)
          p.wet.gibbs[k,g,c(0,cumsum(J))[s]+j,m] <- sum(W.new[[s]][j,month.vec==m]>0,na.rm=TRUE)/sum(!is.na(W.new[[s]][j,month.vec==m]))
          wet.vec <- W.new[[s]][j,month.vec==m & W.new[[s]][j,]>0]
          if (sum(!is.na(wet.vec))>0){
            mean.gibbs[k,g,c(0,cumsum(J))[s]+j,m] <- mean(wet.vec,na.rm=TRUE)
            sd.gibbs[k,g,c(0,cumsum(J))[s]+j,m] <- sd(wet.vec,na.rm=TRUE)
            max.gibbs[k,g,c(0,cumsum(J))[s]+j,m] <- max(wet.vec,na.rm=TRUE)
            
          }
        }
      }
    }
  }
}
post.list <- list(p.wet.gibbs=p.wet.gibbs,mean.gibbs=mean.gibbs,sd.gibbs=sd.gibbs,max.gibbs=max.gibbs)
save(post.list,file=paste(path,"post_list_G5000_042814.RData",sep=""))
load(file=paste(path,"post_list_G5000_042814.RData",sep=""))
for (i in 1:length(post.list)) assign(names(post.list)[i],post.list[[i]])


#save all the data sets for future use, so you don't have to rerun for 5 hours
save.image("~/SkyDrive/IRI/RainfallSimulation/R/Rainfall/post_sim_output042914.RData")


# trace plot function from myRfunctions.R:
tp <- function(gibbs.input,burn=0,end=dim(gibbs.input)[2],nc=3,ylim=range(gibbs.input[,(burn+1):end]),thin=1,...){
  if(nc==1){
    z <- matrix(gibbs.input,ncol=1)
  } else {
    z <- matrix(t(gibbs.input),ncol=nc)
  }
  
  
  G.local <- dim(z)[1]
  #rg <- ifelse(ylim==c(-99,-99),c(min(z[(burn+1):end,]),max(z[(burn+1):end,])),c(ylim[1],ylim[2]))
  #xrg <- ifelse(xlim==c(-99,-99),c(0,length(z[,1])),c(xlim[1],xlim[2]))
  thin.seq <- seq(thin,G.local,by=thin)
  lt <- length(thin.seq)
  plot(thin.seq,z[thin.seq,1],col=1,type="l",ylim=ylim,...,xaxt="n")
  axis(1,at=seq(0,G.local,length=5))
  if(nc > 1){
    for (i in 2:nc) lines(thin.seq,z[thin.seq,i],col=i)
  }
}

############################
#5.Graph Posterior Parameters
############################
#R vs W
pdf(file=paste(path,"fig_R_vs_W.pdf",sep=""),width=7,height=7)
par(mfrow=c(1,1))
for (s in 1:S){
  for (j in 1:J[s]){
    plot(R[[s]][j,R[[s]][j,]>0])
    points(W[[s]][j,is.na(R[[s]][j,])],col=2)
  }
}
dev.off()

# tau
pdf(file=paste(path,"fig_tau.pdf",sep=""),width=7,height=7)
par(mfrow=c(1,1))
tp(tau.gibbs)
#abline(h=true.values$tau,col=4)
dev.off()

# lambda (spatial covariance)
pdf(file=paste(path,"fig_lamba_spatialcov.pdf",sep=""),width=7,height=7)
par(mfrow=c(1,1))
tp(lambda.gibbs)
tp(lambda.gibbs,burn=2500)
#abline(h=true.values$lambda,col=4)
dev.off()

# alpha (t distribution degrees of freedom)
pdf(file=paste(path,"fig_alpha.pdf",sep=""),width=7,height=7)
par(mfrow=c(1,1))
tp(alpha.gibbs)
tp(alpha.gibbs,burn=2500)
#abline(h=true.values$alpha,col=4)
dev.off()

# beta
pdf(file=paste(path,"fig_tobit_trace_beta_P23.pdf",sep=""),width=8,height=6)
par(mfrow=c(2,3))
for (j in 1:P){
  for (s in 1:S){
    tp(beta.gibbs[,,j,s],thin=20,ylim=range(beta.gibbs[,,j,]),main=paste(X.names[j],site.names[s],sep=" "))
    #abline(h=true.values$beta[j,s],col=4,lty=2)
  }
}
dev.off()

# beta.arc
pdf(file=paste(path,"fig_tobit_beta.arc_P23.pdf",sep=""),width=8,height=6)
par(mfrow=c(2,3))
for (s in 1:S){
  tp(beta.arc.gibbs[,,s],main=site.names[s],ylim=range(beta.arc.gibbs),thin=20)
  abline(h=0,lty=2)
  #abline(h=true.values$beta.arc[s],col=4,lwd=2)
}
dev.off()


# Sigma:
pdf(file=paste(path,"fig_tobit_Sigma_P23.pdf",sep=""),width=8,height=6)
for (s in 1:S){
  par(mfrow=c(J[s],J[s]),mar=c(4,3,3,1))
  for (i in 1:J[s]){
    for (j in 1:J[s]){
      tp(Sigma.gibbs[[s]][,,i,j],las=1,thin=20)
      #abline(h=true.values$Sigma[[s]][i,j],col=4,lwd=2)
    }
  }
}
dev.off()

pdf(file=paste(path,"fig_tobit_Sigma_Adiha_P23.pdf",sep=""),width=8,height=6)
par(mfrow=c(1,1),mar=c(5,4,4,2))
s <- 6
for (i in 1:J[s]){
  for (j in i:J[s]){
    tp(Sigma.gibbs[[s]][,,i,j],las=1,thin=20,main=paste(i,j))
  }
}
dev.off()


# mu
pdf(file=paste(path,"fig_tobit_mu_P23.pdf",sep=""),height=8,width=8)
par(mfrow=c(2,2))
for (p in 1:P){
  tp(mu.gibbs[,,p],las=1,main=X.names[p],thin=20)
  #abline(h=true.values$mu[p],col=4)
}
dev.off()


mu.pm <- apply(mu.gibbs[,3001:5000,],3,mean)
z.vec <- X %*% mu.pm

pdf(file=paste(path,"fig_zvector.pdf",sep=""),width=7,height=7)
par(mfrow=c(1,1))
plot(1:(365*2),z.vec[1:(365*2)])
dev.off()

# time trend (January 1st each year for 15 years)
pdf(file=paste(path,"fig_timetrend.pdf",sep=""),width=7,height=7)
plot(seq(1,by=365,length=15),z.vec[seq(1,by=365,length=15)])
dev.off()

# sigma
pdf(file=paste(path,"fig_tobit_trace_sigmalower_P23.pdf",sep=""),height=6,width=9)
par(mfrow=c(2,3))
for (p in 1:P){
  tp(sigma.gibbs[,,p],las=1,thin=20,main=X.names[p])
  #abline(h=true.values$sigma[p],col=4)
}
dev.off()

# mu.arc and tau.arc
pdf(file=paste(path,"fig_tobit_mu.arc.pdf",sep=""),height=6,width=9)
par(mfrow=c(1,2))
tp(mu.arc.gibbs,las=1,thin=20)
#abline(h=true.values$mu.arc,col=4)
dev.off()
pdf(file=paste(path,"fig_tobit_tau.arc.pdf",sep=""),height=6,width=9)
tp(tau.arc.gibbs,las=1,thin=20)
#abline(h=true.values$tau.arc,col=4)
dev.off()


# Look at el Nino effects:
pdf(file=paste(path,"fig_tobit_elnino.pdf",sep=""),height=6,width=9)
se.nino <- numeric(12)
for (i in 1:12) se.nino[i] <- sd(as.numeric(mu.gibbs[,2501:G,i+11]))
par(mfrow=c(1,1))
plot(1:12,mu.pm[12:23],main="Mean monthly El Nino effect",ylim=c(-4,5),xlim=c(0,13))
abline(h=0,lty=2,lwd=2)
text(1:12,mu.pm[12:23],month.names,pos=3)
for (i in 1:12) lines(c(i,i),mu.pm[i+11]+c(-1,1)*se.nino[i])
dev.off()


# Look at mean vector for Z across the whole time span:
pdf(file=paste(path,"fig_tobit_multiyear_mean.pdf",sep=""),width=9,height=7)
par(mfrow=c(1,1))
for (i in 1:6){
  plot(1:(365*3),z.vec[1:(365*3)+(365*3*(i-1))],xaxt="n",xlab="Year",ylab="Mean of Z",type="l",ylim=range(z.vec))
  axis(1,at=seq(0,365*3,length=7),labels=seq(0,by=0.5,length=7)+1992+3*(i-1))
}
dev.off()


# Draw a number of sample paths of multiyear means:
g.seq <- seq(4000,G,length=51)
#ERROR
pdf(file=paste(path,"fig_tobit_multiyear_mean_sample.pdf",sep=""),width=9,height=7)
par(mfrow=c(1,1))
for (i in 1:6){
  temp <- X %*% mu.gibbs[1,g.seq[1],]
  plot(1:(365*3),temp[1:(365*3)+(365*3*(i-1))],xaxt="n",xlab="Year",ylab="Mean of Z",type="l",ylim=c(-50,20),las=1)
  axis(1,at=seq(0,365*3,length=7),labels=seq(0,by=0.5,length=7)+1992+3*(i-1))
  for (g in 2:length(g.seq)){
    temp <- X %*% mu.gibbs[1,g.seq[g],]
    lines(1:(365*3),temp[1:(365*3)+(365*3*(i-1))],col=col.blend)
  }
}
dev.off()

#ERROR
# Now look at the trend without the El Nino effects:
pdf(file=paste(path,"fig_tobit_multiyear_trend.pdf",sep=""),width=9,height=7)
par(mfrow=c(1,1))
for (i in 1:6){
  temp <- X[,1:11] %*% mu.pm[1:11]
  plot(1:(365*3),temp[1:(365*3)+(365*3*(i-1))],xaxt="n",xlab="Year",ylab="Mean of Z",type="l",ylim=c(-50,20),las=1,col=2,lwd=2)
  axis(1,at=seq(0,365*3,length=7),labels=seq(0,by=0.5,length=7)+1992+3*(i-1))
  for (g in 1:length(g.seq)){
    temp <- X[,1:11] %*% mu.gibbs[1,g.seq[g],1:11]
    lines(1:(365*3),temp[1:(365*3)+(365*3*(i-1))],col=col.blend)
  }
  temp <- X[,1:11] %*% mu.pm[1:11]
  lines(1:(365*3),temp[1:(365*3)+(365*3*(i-1))],col=2,lwd=2)
}
dev.off()



#################
#6.Historical Plots
#################
# compute the mean monthly rainfall:
month.vec <- match(months(date.string),month.names)
J.sum <- sum(J)


#KV-Now collect the historical???? data to plot
# From the data:
p.wet.data <- matrix(NA,J.sum,12)
mean.data <- matrix(NA,J.sum,12)
sd.data <- matrix(NA,J.sum,12)
max.data <- matrix(NA,J.sum,12)
for (s in 1:S){
  for (j in 1:J[s]){
    for (m in 1:12){
      #print(c(0,cumsum(J))[s]+j)
      #R is the list of 
      p.wet.data[c(0,cumsum(J))[s]+j,m] <- sum(R[[s]][j,month.vec==m]>0,na.rm=TRUE)/sum(!is.na(R[[s]][j,month.vec==m]))
      wet.vec <- R[[s]][j,month.vec==m & R[[s]][j,]>0]
      mean.data[c(0,cumsum(J))[s]+j,m] <- mean(wet.vec,na.rm=TRUE)
      sd.data[c(0,cumsum(J))[s]+j,m] <- sd(wet.vec,na.rm=TRUE)
      max.data[c(0,cumsum(J))[s]+j,m] <- max(wet.vec,na.rm=TRUE)
      #print(paste(s,j,m,length(wet.vec),sum(!is.na(wet.vec))),sep=",")
    }
  }
}

#ERROR
max.data[max.data==-Inf] <- NA

#KV-Plot the actual original data:
par(mfrow=c(1,1))


# probability of rain:
quartz()
pdf(file=paste(path,"PercWetDays_Historical.pdf",sep=""),height=6,width=9)
plot(-10,-10,xlim=c(1,12),ylim=c(0,1),las=1,xlab="Month",ylab="% Wet Days",xaxt="n",main="% of wet days, HISTORICAL")
axis(1,at=1:12,labels=substr(month.names,1,1))
for (i in 1:15) lines(1:12,p.wet.data[i,],col=i)
dev.off()


# mean wet day rainfall:

pdf(file=paste(path,"MeanWetDays_Historical.pdf",sep=""),height=6,width=9)
plot(-10,-10,xlim=c(1,12),ylim=c(0,max(mean.data,na.rm=TRUE)),las=1,xlab="Month",ylab="Mean (mm)",xaxt="n",mean="Mean wet day rainfall,HISTORICAL")
axis(1,at=1:12,labels=substr(month.names,1,1))
for (i in 1:15) lines(1:12,mean.data[i,],col=i)
dev.off()


# sd wet day rainfall:
quartz()
pdf(file=paste(path,"SdWetDays_Historical.pdf",sep=""),height=6,width=9)
plot(-10,-10,xlim=c(1,12),ylim=c(0,max(sd.data,na.rm=TRUE)),las=1,xlab="Month",ylab="Sd (mm)",xaxt="n",main="Sd of wet day rainfall,HISTORICAL")
axis(1,at=1:12,labels=substr(month.names,1,1))
for (i in 1:15) lines(1:12,sd.data[i,],col=i)
dev.off()


# max wet day rainfall:
pdf(file=paste(path,"MaxWetDays_Historical.pdf",sep=""),height=6,width=9)
plot(-10,-10,xlim=c(1,12),ylim=c(0,max(max.data,na.rm=TRUE)),las=1,xlab="Month",ylab="Max (mm)",xaxt="n",main="Max of wet day rainfall,HISTORICAL")
axis(1,at=1:12,labels=substr(month.names,1,1))
for (i in 1:15) lines(1:12,max.data[i,],col=i)


#################
#7.Posterior Plots
#################
# Plot  posterior predictive checks:
#use burn in data
pb <- seq(3000,5000,by=5)
pdf(file=paste(path,"fig_tobit_ppred_max.pdf",sep=""),width=10,height=8)
par(mfrow=c(3,4))
for (i in 1:J.sum){
  for (m in 1:12){
    #g.out <- p.wet.gibbs[,pb,i,m]; d.out <- p.wet.data[i,m]
    #lab <- "P(Wet)"
    #g.out <- mean.gibbs[,pb,i,m]; d.out <- mean.data[i,m]
    #lab <- "Mean wet days (mm)"
    #g.out <- sd.gibbs[,pb,i,m]; d.out <- sd.data[i,m]
    #lab <- "Sd wet days (mm)"
    g.out <- max.gibbs[,pb,i,m]; d.out <- max.data[i,m]
    lab <- "Max rainfall (mm)"
    if (!is.na(d.out)){
      xrg <- range(c(g.out,d.out),na.rm=TRUE)
      hist(g.out,main=month.names[m],las=1,xlab=lab,xlim=xrg)
      abline(v=d.out,col=2,lwd=2)
    } else {
      plot(-10,-10,xlim=c(0,1),ylim=c(0,1))
      text(0.5,0.5,"No Data")
    }
  }
}
dev.off()


#KV-More posterior predictive checks:
# Compute p-values for each of the 4 x 15 x 12 checks:
p.val <- array(NA,dim=c(4,J.sum,12))
for (i in 1:J.sum){
  for (m in 1:12){
    p.val[1,i,m] <- sum(p.wet.data[i,m] > p.wet.gibbs[,pb,i,m])/(K*length(pb))
    p.val[2,i,m] <- sum(mean.data[i,m] > mean.gibbs[,pb,i,m])/(K*length(pb))
    p.val[3,i,m] <- sum(sd.data[i,m] > sd.gibbs[,pb,i,m])/(K*length(pb))
    p.val[4,i,m] <- sum(max.data[i,m] > max.gibbs[,pb,i,m])/(K*length(pb))
  }
}

#KV plot comparable posterior graphs to historical
post <- p.wet.gibbs[1,4001:5000,1:15,1:12]
post_p <- colMeans(post)

pdf(file=paste(path,"PercWetDays_Posterior.pdf",sep=""),height=6,width=9)
plot(-10,-10,xlim=c(1,12),ylim=c(0,1),las=1,xlab="Month",ylab="% Wet Days",xaxt="n",main="% of wet days, HISTORICAL")
axis(1,at=1:12,labels=substr(month.names,1,1))
for (i in 1:15) lines(1:12,post_p[i,],col=i)
dev.off()

# mean wet day rainfall:
post <- mean.gibbs[1,4001:5000,1:15,1:12]
post_m <- colMeans(post)

pdf(file=paste(path,"MeanWetDays_Posterior.pdf",sep=""),height=6,width=9)
plot(-10,-10,xlim=c(1,12),ylim=c(0,max(post_m,na.rm=TRUE)),las=1,xlab="Month",ylab="Mean (mm)",xaxt="n",mean="Mean wet day rainfall,HISTORICAL")
axis(1,at=1:12,labels=substr(month.names,1,1))
for (i in 1:15) lines(1:12,post_m[i,],col=i)
dev.off()


# sd wet day rainfall:
post <- sd.gibbs[1,4001:5000,1:15,1:12]
post_s <- colMeans(post)

pdf(file=paste(path,"SdWetDays_Posterior.pdf",sep=""),height=6,width=9)
plot(-10,-10,xlim=c(1,12),ylim=c(0,max(post_s,na.rm=TRUE)),las=1,xlab="Month",ylab="Sd (mm)",xaxt="n",main="Sd of wet day rainfall,HISTORICAL")
axis(1,at=1:12,labels=substr(month.names,1,1))
for (i in 1:15) lines(1:12,post_s[i,],col=i)
dev.off()


# max wet day rainfall:
post <- max.gibbs[1,4001:5000,1:15,1:12]
post_ma <- colMeans(post)

pdf(file=paste(path,"MaxWetDays_Posterior.pdf",sep=""),height=6,width=9)
plot(-10,-10,xlim=c(1,12),ylim=c(0,max(post_ma,na.rm=TRUE)),las=1,xlab="Month",ylab="Max (mm)",xaxt="n",main="Max of wet day rainfall,HISTORICAL")
axis(1,at=1:12,labels=substr(month.names,1,1))
for (i in 1:15) lines(1:12,post_ma[i,],col=i)
dev.off()









#KV-Plot the distribution of the pvalues for all the moments from the posterior-simulated W's
pdf(file=paste(path,"fig_pvals_hist.pdf",sep=""),width=10,height=8)
stat.names <- c("P(wet)","Mean","Sd","Max")
par(mfrow=c(2,2))
for (i in 1:4) hist(p.val[i,,],xlim=c(0,1),breaks=seq(0,1,0.1),main=stat.names[i])
dev.off()


# Plot some values of W against the observed values (R):
pdf(file=paste(path,"fig_ppred_WvsObservedR.pdf",sep=""),width=10,height=8)
par(mfrow=c(2,1))
plot(W.new[[1]][1,na.mat[[1]][1,]==0])
plot(R[[1]][1,na.mat[[1]][1,]==0],col=2)
dev.off()



# trace plot of covariance matrix elements:
pdf(file=paste(path,"fig_ppred_cov_matrix.pdf",sep=""),width=10,height=8)
s <- 1
par(mfrow=c(2,2))
tp(Sigma.gibbs[[s]][,,1,1],las=1)
tp(Sigma.gibbs[[s]][,,1,2],las=1)
tp(Sigma.gibbs[[s]][,,2,1],las=1)
tp(Sigma.gibbs[[s]][,,2,2],las=1)
dev.off()


# Compute some MLEs:
#ERROR no true.values
Sigma.hat <- as.list(rep(NA,S))
for (s in 1:S){
  mu.W <- t(matrix(true.values$Z[,s],T,J[s])) - matrix(true.values$beta.arc[s]*X.arc[[s]],J[s],T)
  Sigma.hat[[s]] <- var(t(true.values$W[[s]] - mu.W))
}


ba.mle <- numeric(S)
n.ba <- numeric(S)
for (s in 1:S){
  #ba.mle[s] <- mean(true.values$W[[s]][1,] - true.values$W[[s]][2,])
  sel <- R[[s]][1,] > 0 & R[[s]][2,] > 0
  ba.mle[s] <- mean(R[[s]][1,sel] - R[[s]][2,sel],na.rm=TRUE)
  n.ba[s] <- sum(sel,na.rm=TRUE)
}


#s <- 1
pdf(file=paste(path,"fig_w_samp",s,".pdf",sep=""),width=9,height=5)
par(mfrow=c(1,2))
for (i in 1:n.samp){
  for (j in 1:J[s]){
    tp(W.gibbs[[s]][,1:ng,j,i],main=paste(date.string[w.samp[i]%%365],round(R[[s]][j,w.samp[i]],2),sep=": "),ylim=range(W.gibbs[[s]][,,,i]))
    abline(h=0,lty=2)
  }
}
dev.off()



pdf(file=paste(path,"fig_z_samp.pdf",sep=""),width=9,height=7)
par(mfrow=c(2,2))
for (i in 1:n.samp){
  plot(1:ng,Z.gibbs[1,1:ng,i,1],type="l",ylim=range(Z.gibbs[,1:ng,i,1]))
  for (k in 2:3) lines(1:ng,Z.gibbs[k,1:ng,i,1],col=k)
  abline(h=true.values$Z[w.samp[i],1])
}
dev.off()


# Try adjustments to mu and beta (later...)
m.gibbs <- array(NA,dim=c(K,G,P))
for (k in 1:K){
  for (p in 1:P){
    m.gibbs[k,,p] <- mean(beta.gibbs[k,,p,])
  }
}

b.gibbs <- array(NA,dim=c(K,G,P,S))
for (p in 1:P){
  for (s in 1:S){
    b.gibbs[,,p,s] <- beta.gibbs[,,p,s] - mean(beta.gibbs[,,p,s])
  }
}
