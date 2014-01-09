###############################################################################################################
# Read in the data
rm(list=ls())
#source("~/Stats/Misc/myRfunctions.R")
#data.path <- "~/Stats/IndexInsurance/AdiHa supporting data/"
setwd("/Users/kathrynvasilaky/SkyDrive/IRI/RainfallSimulation/R")
#setwd("/Users/katyav/SkyDrive/IRI/RainfallSimulation/R")
path <- getwd()#"~/SkyDrive/IRI/RainfallSimulation/R/"  # enter in here wherever you want to store the scripts and data files
path <- paste(path,'/', sep='')
source(paste(path,"R code multisite covariance scripts.R",sep=""))  # read in some scripts I wrote
#library(bayesm)
library(MASS)

# load in old data set with Adi Ha short data sets and a couple others:
load(paste(path,"ethiopia_full_data.RData",sep=""))  # read in the data, saved as an R object
#KV-15 rows or data from 15 different sites, reorders them in 5's?
data <- data[c(6,1,7,2,8,3,9,4,10,5,11,12,13,14,15),]
T <- dim(data)[2]
site.mat <- cbind(c(1,1,2,2,3,3,4,4,5,5,6,6,6,6,6),c(1,2,1,2,1,2,1,2,1,2,1,2,3,4,5))
arc <- c(1,0,1,0,1,0,1,0,1,0,1,0,0,0,0)

S <- 6
L <- c(2,2,2,2,2,5)
L.sum <- sum(L)
month.days <- c(31,28,31,30,31,30,31,31,30,31,30,31)
month <- c(rep(rep(1:12,month.days),49),rep(1:12,month.days)[1:209])

# 1992 onward:
data <- data[,(365*31+1):T] # 1992 onward, all 15 sites:
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
nino <- matrix(t(as.matrix(nino[,-1])),ncol=1)
nino.dates <- paste(rep(1871:2010,each=12),unique(months(date.string)),sep=" ")
nino <- nino[1:1672]; nino.dates <- nino.dates[1:1672]  # nino only goes through April 2010.

# subselect just the months that align with the rainfall data:
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

# Enter ARC indicator information:
# KV-ARC data is http://journals.ametsoc.org/doi/abs/10.1175/JCLI-D-12-00206.1
# KV-"(ARC) project aims to create an independent climate data record of sea surface 
#temperatures (SSTs) covering recent decades that can be used for climate change analysis"
X.arc <- as.list(rep(NA,S))
for (s in 1:S){
  X.arc[[s]] <- matrix(0,J[s],1)
  X.arc[[s]][1,1] <- 1
}

# Real data:
# KV - A list of lists ... the same dimension as na.mat
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

# Set starting values for the MCMC
K <- 3


######## skip this if using real data ###########

######## resume here if using real data ###########

### Wei, you can stop here -- the rest is setting up the MCMC, and I have to give you two C functions, "zdraw" and "draw_gamma", to do the MCMC.
### Hopefully the above will at least give you a start on learning how to fit a model like this one.

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



### Priors
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



# Set up storage for gibbs samples:
K <- 3
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

# Start the MCMC:
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
save(gibbs.list,file=paste(path,"gibbs_out_11112013_G5000.RData",sep=""))


load(file=paste(path,"gibbs_out_11112013_G5000.RData",sep=""))
for (i in 1:length(gibbs.list)) assign(names(gibbs.list)[i],gibbs.list[[i]])


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
pdf(file=paste(path,"fig_lamba.pdf",sep=""),width=7,height=7)
par(mfrow=c(1,1))
tp(lambda.gibbs)
tp(lambda.gibbs,burn=2500)
#abline(h=true.values$lambda,col=4)
dev.off()

# alpha (t distribution degrees of freedom)
pdf(file=paste(path,"fig_alpha.pdf",sep=""),width=7,height=7)
par(mfrow=c(1,1))
tp(alpha.gibbs)
tp(alpha.gibbs,burn=1000)
abline(h=true.values$alpha,col=4)
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








# New starting points:
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
save(start.list,file=paste(path,"start_list_v3.RData",sep=""))
# v3 is the G=5000 output from 3/12 that includes P=23 predictors and alpha=10








# Posterior predictive checks:

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

# test case:
k <- 1
g <- 75

# Collect Sigma into a list:
Sigma.sim <- as.list(rep(NA,S))
for (s in 1:S) Sigma.sim[[s]] <- Sigma.gibbs[[s]][k,g,,]

# Simulate rainfall:
W.new <- sim.W(alpha=alpha.gibbs[k,g],beta=beta.gibbs[k,g,,],lambda=lambda.gibbs[k,g],tau=tau.gibbs[k,g],beta.arc=beta.arc.gibbs[k,g,],
               Sigma=Sigma.sim,X=X,X.arc=X.arc,na.preserve=TRUE,na.mat=na.mat)

# compute the mean monthly rainfall:
month.vec <- match(months(date.string),month.names)
J.sum <- sum(J)

# From the data:
p.wet.data <- matrix(NA,J.sum,12)
mean.data <- matrix(NA,J.sum,12)
sd.data <- matrix(NA,J.sum,12)
max.data <- matrix(NA,J.sum,12)
for (s in 1:S){
  for (j in 1:J[s]){
    for (m in 1:12){
      #print(c(0,cumsum(J))[s]+j)
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

# Plot the data:
par(mfrow=c(1,1))

# probability of rain:
plot(-10,-10,xlim=c(1,12),ylim=c(0,1),las=1,xlab="Month",ylab="% Wet Days",xaxt="n",main="% of wet days")
axis(1,at=1:12,labels=substr(month.names,1,1))
for (i in 1:15) lines(1:12,p.wet.data[i,],col=i)

# mean wet day rainfall:
plot(-10,-10,xlim=c(1,12),ylim=c(0,max(mean.data,na.rm=TRUE)),las=1,xlab="Month",ylab="Mean (mm)",xaxt="n",mean="Mean wet day rainfall")
axis(1,at=1:12,labels=substr(month.names,1,1))
for (i in 1:15) lines(1:12,mean.data[i,],col=i)

# sd wet day rainfall:
plot(-10,-10,xlim=c(1,12),ylim=c(0,max(sd.data,na.rm=TRUE)),las=1,xlab="Month",ylab="Sd (mm)",xaxt="n",main="Sd of wet day rainfall")
axis(1,at=1:12,labels=substr(month.names,1,1))
for (i in 1:15) lines(1:12,sd.data[i,],col=i)

# max wet day rainfall:
plot(-10,-10,xlim=c(1,12),ylim=c(0,max(max.data,na.rm=TRUE)),las=1,xlab="Month",ylab="Max (mm)",xaxt="n",main="Max of wet day rainfall")
axis(1,at=1:12,labels=substr(month.names,1,1))
for (i in 1:15) lines(1:12,max.data[i,],col=i)





# From the gibbs output:
p.wet.gibbs <- array(NA,dim=c(K,G,J.sum,12)) # monthly percentage of wet days:
mean.gibbs <- array(NA,dim=c(K,G,J.sum,12)) # monthly percentage of wet days:
sd.gibbs <- array(NA,dim=c(K,G,J.sum,12)) # monthly percentage of wet days:
max.gibbs <- array(NA,dim=c(K,G,J.sum,12)) # monthly percentage of wet days:
for (k in 1:K){
  for (g in 1:G){
    if (g%%20==0) print(g)
    # Collect Sigma into a list:
    Sigma.sim <- as.list(rep(NA,S))
    for (s in 1:S) Sigma.sim[[s]] <- Sigma.gibbs[[s]][k,g,,]
    # Simulate rainfall:
    W.new <- sim.W(alpha=alpha.gibbs[k,g],beta=beta.gibbs[k,g,,],lambda=lambda.gibbs[k,g],tau=tau.gibbs[k,g],beta.arc=beta.arc.gibbs[k,g,],
                   Sigma=Sigma.sim,X=X,X.arc=X.arc,na.preserve=TRUE,na.mat=na.mat)
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
save(post.list,file=paste(path,"post_list_G5000_0312.RData",sep=""))
load(file=paste(path,"post_list_G5000_0312.RData",sep=""))
for (i in 1:length(post.list)) assign(names(post.list)[i],post.list[[i]])



# Plot these posterior predictive checks:
pb <- seq(3000,5000,by=5)
pdf(file=paste(path,"fig_tobit_ppred.pdf",sep=""),width=10,height=8)
par(mfrow=c(3,4))
for (i in 1:J.sum){
  for (m in 1:12){
    g.out <- p.wet.gibbs[,pb,i,m]; d.out <- p.wet.data[i,m]
    lab <- "P(Wet)"
    #g.out <- mean.gibbs[,pb,i,m]; d.out <- mean.data[i,m]
    #lab <- "Mean wed days (mm)"
    #g.out <- sd.gibbs[,pb,i,m]; d.out <- sd.data[i,m]
    #lab <- "Sd wet days (mm)"
    #g.out <- max.gibbs[,pb,i,m]; d.out <- max.data[i,m]
    #lab <- "Max rainfall (mm)"
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


# Compute p-values for each of the 4 x 15 x 12 checks:
#pdf(file=paste(path,"fig_pvals.pdf",sep=""),width=10,height=8)
p.val <- array(NA,dim=c(4,J.sum,12))
for (i in 1:J.sum){
  for (m in 1:12){
    p.val[1,i,m] <- sum(p.wet.data[i,m] > p.wet.gibbs[,pb,i,m])/(K*length(pb))
    p.val[2,i,m] <- sum(mean.data[i,m] > mean.gibbs[,pb,i,m])/(K*length(pb))
    p.val[3,i,m] <- sum(sd.data[i,m] > sd.gibbs[,pb,i,m])/(K*length(pb))
    p.val[4,i,m] <- sum(max.data[i,m] > max.gibbs[,pb,i,m])/(K*length(pb))
  }
}
#dev.off()

pdf(file=paste(path,"fig_pvals_hist.pdf",sep=""),width=10,height=8)
stat.names <- c("P(wet)","Mean","Sd","Max")
par(mfrow=c(2,2))
for (i in 1:4) hist(p.val[i,,],xlim=c(0,1),breaks=seq(0,1,0.1),main=stat.names[i])
dev.off()

# Plot these posterior predictive checks:
pb <- seq(3000,5000,by=5)
pdf(file=paste(path,"fig_tobit_ppred.pdf",sep=""),width=10,height=8)
par(mfrow=c(3,4))
for (i in 1:J.sum){
  for (m in 1:12){
    g.out <- p.wet.gibbs[,pb,i,m]; d.out <- p.wet.data[i,m]
    lab <- "P(Wet)"
    g.out <- mean.gibbs[,pb,i,m]; d.out <- mean.data[i,m]
    lab <- "Mean wed days (mm)"
    g.out <- sd.gibbs[,pb,i,m]; d.out <- sd.data[i,m]
    lab <- "Sd wet days (mm)"
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


s <- 1
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

















##### Old junk below here...

# EOF

# load bayesm library
#set.seed(66)
#rmvst(nu=5,mu=c(rep(0,2)),root=chol(matrix(c(2,1,1,2),ncol=2)))

n.sim <- 5000
df <- 5

x1 <- matrix(NA,n.sim,2)
x2 <- mvrnorm(n.sim,mu=c(0,0),Sigma=matrix(c(1,0.8,0.8,1),2,2))
x3 <- matrix(NA,n.sim,2)
gamma.vec <- rgamma(n.sim,shape=df/2,scale=2/df)
x4 <- mvrnorm(n.sim,mu=c(0,0),Sigma=matrix(c(1,0.8,0.8,1),2,2))/sqrt(gamma.vec)
for (i in 1:n.sim){
  x1[i,] <- rmvst(nu=df,mu=c(rep(0,2)),root=chol(matrix(c(1,0.8,0.8,1),ncol=2)))
  gamma <- rgamma(1,shape=df/2,scale=2/df)
  x3[i,] <- mvrnorm(n=1,mu=rep(0,2),Sigma=matrix(c(1,0.8,0.8,1),2,2)/gamma)
}
par(mfrow=c(2,2))
rg <- c(-15,15)
plot(x2[,1],x2[,2],xlim=rg,ylim=rg)
plot(x1[,1],x1[,2],xlim=rg,ylim=rg)
plot(x3[,1],x3[,2],xlim=rg,ylim=rg)
plot(x4[,1],x4[,2],xlim=rg,ylim=rg)





# Draw random multivariate student-t variables:
W <- as.list(rep(NA,S))
for (s in 1:S){
  print(s)
  W[[s]] <- matrix(NA,J[s],T)
  X.w <- cbind(rep(1,J[s]),X.arc[[s]])
  for (t in 1:T){
    W[[s]][,t] <- rmvst(nu=alpha,mu=X.w %*% rbind(Z[s,t],beta.arc[s]),root=chol(Sigma[[s]]))
  }
}

# Plot the mean signal:
#x.mat <- cbind(rep(1,365*3),sin.mat[1:(365*3),1:2],cos.mat[1:(365*3),1:2])
#plot(1:(365*3), x.mat %*% matrix(beta,ncol=1))
# Good, it looks how it should.

w.mis <- as.list(rep(NA,S))
n.mis <- as.list(rep(NA,S))
w.dry <- as.list(rep(NA,S))
n.dry <- as.list(rep(NA,S))
for (s in 1:S){
  w.mis[[s]] <- matrix(NA,J[s],T)
  w.dry[[s]] <- matrix(NA,J[s],T)
  n.mis[[s]] <- numeric(J[s])
  n.dry[[s]] <- numeric(J[s])
  for (j in 1:J[s]){
    w.mis[[s]][j,] <- is.na(R[[s]][j,])
    w.dry[[s]][j,] <- R[[s]][j,]==0
    n.mis[[s]][j] <- sum(is.na(R[[s]][j,]))
    n.dry[[s]][j] <- sum(R[[s]][j,]==0,na.rm=TRUE)
  }
}

# Test a few distributions:
n <- 50000
mu <- 10
sigma <- 5
alpha <- 10
gamma <- rgamma(n,shape=alpha/2,rate=alpha/2)
y <- rnorm(n,mu,sqrt(sigma^2/gamma))

plot(density(y))
lines(seq(-20,40,0.01),dst(seq(-20,40,0.01),df=alpha,mu=10,sigma=5),col=2)
# Good, the mixture works.


# Sort out gamma vs. inverse gamma:
# density of an inverse gamma with shape and scale parameters
dig <- function(x,shape,scale) scale^shape/gamma(shape)*x^(-(shape+1))*exp(-scale/x)


y <- rgamma(n,shape=5,scale=2)
plot(density(1/y))
lines(seq(0,2,0.01),dig(seq(0,2,0.01),shape=5,scale=1/2),col=2)
# OK: if x ~ gamma(shape=alpha,scale=beta), then 1/x ~ inverse-gamma(shape=alpha,scale=1/beta)




# Old attempt at sampling gamma:
t1 <- Sys.time()
for (s in 1:S){      
  mu.W <- t(Z[,s] + matrix(rep(X.arc[[s]]*beta.arc[s],each=T),ncol=J[s]))
  #prac <- t(W[[s]]-mu.W) %*% (W[[s]]-mu.W)
  #W[[s]][,1] %*% matrix(W[[s]][,1],nrow=1)*solve(Sigma[[s]])
  gamma.mat[,s] <- rgamma(T,shape=(J[s]+alpha)/2,rate=(diag(t(W[[s]]-mu.W) %*% solve(Sigma[[s]]) %*% (W[[s]]-mu.W))+alpha)/2)
}
t2 <- Sys.time()
t2-t1


# Old function to draw Sigma:
draw.Sigma <- function(W,Z,X.arc,beta.arc,gamma.mat){
  Sigma.draw <- as.list(rep(NA,S))
  for (s in 1:S){
    mu.mat <- t(matrix(Z[,s],T,J[s]) + matrix(rep(X.arc[[s]]*beta.arc[s],each=T),ncol=J[s]))  # J[s] x T matrix
    Sigma.draw[[s]] <- rwishart(nu=dim(W[[s]])[2] + J[s],V=solve(diag(J[s]) + (W[[s]]-mu.mat)%*%t(W[[s]]-mu.mat)))$IW
  }
  return(Sigma.draw)
}

ch.R <- t(chol(prior.inv))
dyn.load(paste(path,"z_draw.so",sep=""))
is.loaded("z_draw")
ch <- .C("z_draw",a_in=as.double(prior.inv),n_in=as.integer(S),p=as.double(rep(0,S)))
ch.C <- matrix(0,S,S)
ch.C[lower.tri(ch.C)] <- matrix(ch$a_in,S,S)[lower.tri(ch.C)]
diag(ch.C) <- ch$p
ch.C
# good, it works!



# Draw Z | W, gamma, beta.arc, Sigma, beta, lambda, tau:
t1 <- Sys.time()
prior.cov <- solve(R.cov(lambda,d.mat))/tau^2
sigma.z <- numeric(s)
for (s in 1:S) sigma.z[s] <- 1/sum(solve(Sigma[[s]]))
Sig.z.mat <- 1/t(gamma.mat)*sigma.z
mu.z <- matrix(0,S,T)
pc.xbeta <- prior.cov %*% t(X%*%beta)
for (s in 1:S){
  mu.z[s,] <- rep(-1,J[s]) %*% solve(Sigma[[s]]) %*% (matrix(X.arc[[s]]*beta.arc[s],J[s],T) - W[[s]])
  for (t in 1:T){
    sz.inv <- solve(diag(Sig.z.mat[,t]))
    A.inv <- solve(sz.inv + prior.cov)
    mu.t <- A.inv %*% (sz.inv %*% mu.z[,t] + pc.xbeta[,t])
    Z[t,] <- mvrnorm(n=1,mu=mu.t,Sigma=A.inv)
  }
}
t2 <- Sys.time()
t2-t1


for (t in 1:T){
  rate.vec[t] <- (matrix(W[[s]][,t]-mu.W[,t],nrow=1) %*% S.inv %*% matrix(W[[s]][,t]-mu.W[,t],ncol=1) + alpha)/2
}


# Draw Z in R (without the C function call) to make sure its working right.
if (FALSE){
  #t1 <- Sys.time()
  sigma.z <- numeric(s)
  for (s in 1:S) sigma.z[s] <- sum(solve(Sigma[[s]]))
  mu.z <- matrix(0,S,T)
  for (s in 1:S) mu.z[s,] <- (matrix(-1,1,J[s]) %*% solve(Sigma[[s]]) %*% (matrix(X.arc[[s]]*beta.arc[s],J[s],T) - W[[s]]))/sigma.z[s]
  prior.inv <- solve(R.cov(lambda,d.mat))/tau^2
  Pxb <- prior.inv %*% t(X%*%beta)
  Sig <- t(gamma.mat)*sigma.z
  # Store the mean and covariance matrix:
  #A.inv <- array(NA,dim=c(T,S,S))
  #mu.all <- matrix(NA,S,T)
  for (t in 1:T){
    if (t%%2000==0) print(t)
    sz.inv <- diag(Sig[,t])
    #A.inv[t,,] <- solve(sz.inv + prior.cov)
    #mu.all[,t] <- A.inv[t,,] %*% (sz.inv %*% mu.z[,t] + pc.xbeta[,t])
    A.inv <- solve(sz.inv + prior.inv) # A.inv is the covariance of the new draw of z.
    mu.all <- A.inv %*% (sz.inv %*% mu.z[,t] + Pxb[,t])
    Z[t,] <- mvrnorm(n=1,mu=mu.all,Sigma=A.inv)
  }
  #t2 <- Sys.time()
  #t2-t1
}


prac <- matrix(NA,1000,6)
for (i in 1:1000){
  if (i%%200==0) print(i)
  prac[i,] <- t(matrix(z.draw$Smu_in,S,T))[1,]
}





# To do:
# (1) check generation of W in the simulated data.
# Done, and was correct
# (2) check the variance of Z in the Gibbs draw of Z
# (3) check the mean of Z in the Gibbs draw
# Done, and fixed
# (4) Check MLE for correlations and beta.arc from data to see if we're getting it right after all
