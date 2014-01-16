#KV-This is a hack job at creating code to simulate observations from Aug 2010 to Aug 2011, using the save Gibss samples. 
#KV-You need to run R_code_tobit_real1.R, or upload the saves parameters below from the last run completed. 

#load parameter samples
load(file=paste(path,"gibbs_out_01152014_G5000.RData",sep=""))
for (i in 1:length(gibbs.list)) assign(names(gibbs.list)[i],gibbs.list[[i]])


##################
#Create new X, seasonal data for Aug 1,2010 until July 31, 2011
#I at August 1st, else it's kind of awkard w/ 2 days in July
T<-6781
num_days_to_sim <- 365
day_range <- (T+1):(T+num_days_to_sim)
# Create periodic predictor variables:
# KV - for each column (10 all together) he is creating a sine 
# and a cosine
# In one year, the wave goes up and down once, plot(sin.mat[,1])
M <- 10
m.seq <- c(1:10)
sin.mat <- matrix(NA,num_days_to_sim,M)
cos.mat <- matrix(NA,num_days_to_sim,M)
for (m in 1:M){
  sin.mat[,m] <- sin(2*pi*(day_range)*m.seq[m]/365)
  cos.mat[,m] <- cos(2*pi*(day_range)*m.seq[m]/365)
}


#################################
#THIS CHUNK USES NEW EL NINO DATA
#AND RECONSTRUCTS THE NINO MATRIX
#################################
# Read in Nino 3.4 index:
# KV - this index may come from here: http://www.cgd.ucar.edu/cas/catalog/climind/TNI_N34/
#nino <- read.table(paste(path,"nino34.long.data",sep=""),colClasses=rep("numeric",13))
nino <- read.table(paste(path,"nin34_long_data_new.data",sep=""),colClasses=rep("numeric",13))
nino <- matrix(t(as.matrix(nino[,-1])),ncol=1)
#Old date string
date.string <- as.Date(1:(T+5)-1,origin="1992-01-01")
leap.year.index <- 60 + c(0,365*4+1,365*8+2,365*12+3,365*16+4)
date.string <- date.string[-leap.year.index]
month.names <- unique(months(date.string))

nino.dates <- paste(rep(1871:2011,each=12),unique(months(date.string)),sep=" ")
#nino <- nino[1:1672]; nino.dates <- nino.dates[1:1672]  # nino only goes through April 2010.
# subselect just the months that align with the rainfall data:
w <- which(nino.dates=="2010 August") # this is where the real data end
e<- which(nino.dates=="2011 July") # this is where the real data end
nino <- nino[w:e]
nino.dates <- nino.dates[w:e]
# OK, now nino starts in August 2010 and ends August 2011
# center nino at its mean
nino <- nino - mean(nino)

#create NEW DATE STRING starting Aug 2010
date.string <- as.Date(1:365,origin="2010-07-31")
month.names <- unique(months(date.string))
#Biz as usual
month.vec <- match(months(date.string),month.names)
month.mat <- matrix(0,num_days_to_sim,12)
#month.mat[cbind(1:T,month.vec)] <- rep(nino[1:223],c(rep(month.days,18),month.days[1:6],28))
month.mat[cbind(1:num_days_to_sim,month.vec)] <- rep(nino[1:12],month.days)

#############################################

#use 
#create fake el nino matrix
#month.mat.fake<-matrix(.5,365,12)

# Create NEW X_t, which is the same across sites:
#X <- cbind(rep(1,num_days_to_sim), (day_range)/T-0.5, ((day_range)/T-0.5)^2,sin.mat[,1:4],cos.mat[,1:4],month.mat)
X11 <- cbind(rep(1,num_days_to_sim), (day_range+T+1)/T-0.5, ((day_range+T+1)/T-0.5)^2,sin.mat[,1:4],cos.mat[,1:4],month.mat.fake)

#test is 365x12
#X is 6779x23
#beta is 23x6
#beta.arc is 1x6

# changed the names here from sim.W to sim11.W etc.
# function to simulate rainfall given parameters from gibbs output:

sim.W <- function(alpha,beta,lambda,tau,beta.arc,Sigma,X11,X.arc,na.preserve=TRUE,na.mat){
  Z.sim <- mvrnorm(T,rep(0,S),tau^2*R.cov(lambda,d.mat)) + X11 %*% beta
  gamma.mat <- matrix(rgamma(T*S,shape=alpha/2,scale=2/alpha),T,S)  
  W.new <- as.list(rep(NA,S))
  for (s in 1:S){
    mu.W <- t(Z.sim[,s] + t(matrix(X.arc[[s]]*beta.arc[s],J[s],T)))
    temp <- t( mvrnorm(n=T,mu=rep(0,J[s]),Sigma=Sigma[[s]])/sqrt(gamma.mat[,s]) )
    temp <- temp + mu.W
    temp[temp<0] <- 0
    W.new[[s]] <- temp
    #if (na.preserve) W.new[[s]][na.mat[[s]]==1] <- NA
  }
  #return(unlist(W.new)[unlist(na.mat)==0])
  return(W.new)
}


#I want to only sample from the end of the chain



# test case:
k <- 1
#use the end of the chain to simulate data
g <- 5000

# Collect Sigma into a list:
Sigma.sim <- as.list(rep(NA,S))
for (s in 1:S) Sigma.sim[[s]] <- Sigma.gibbs[[s]][k,g,,]

# Simulate rainfall:
W.new <- sim.W(alpha=alpha.gibbs[k,g],beta=beta.gibbs[k,g,,],lambda=lambda.gibbs[k,g],tau=tau.gibbs[k,g],beta.arc=beta.arc.gibbs[k,g,],
               Sigma=Sigma.sim,X=X11,X.arc=X.arc,na.preserve=TRUE,na.mat=na.mat)

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
pdf(file=paste(path,"Predicted2011ProbRainfall.pdf",sep=""),width=10,height=8)
plot(-10,-10,xlim=c(1,12),ylim=c(0,1),las=1,xlab="Month",ylab="% Wet Days",xaxt="n",main="% of wet days")
axis(1,at=1:12,labels=substr(month.names,1,1))
for (i in 1:15) lines(1:12,p.wet.data[i,],col=i)
dev.off()

# mean wet day rainfall:
pdf(file=paste(path,"Predicted2011MeanRainfall.pdf",sep=""),width=10,height=8)
plot(-10,-10,xlim=c(1,12),ylim=c(0,max(mean.data,na.rm=TRUE)),las=1,xlab="Month",ylab="Mean (mm)",xaxt="n",main="Mean wet day rainfall")
axis(1,at=1:12,labels=substr(month.names,1,1))
for (i in 1:15) lines(1:12,mean.data[i,],col=i)
dev.off()

# sd wet day rainfall:
pdf(file=paste(path,"Predicted2011SDRainfall.pdf",sep=""),width=10,height=8)
plot(-10,-10,xlim=c(1,12),ylim=c(0,max(sd.data,na.rm=TRUE)),las=1,xlab="Month",ylab="Sd (mm)",xaxt="n",main="Sd of wet day rainfall")
axis(1,at=1:12,labels=substr(month.names,1,1))
for (i in 1:15) lines(1:12,sd.data[i,],col=i)
dev.off()

# max wet day rainfall:
pdf(file=paste(path,"Predicted2011MaxRainfall.pdf",sep=""),width=10,height=8)
plot(-10,-10,xlim=c(1,12),ylim=c(0,max(max.data,na.rm=TRUE)),las=1,xlab="Month",ylab="Max (mm)",xaxt="n",main="Max of wet day rainfall")
axis(1,at=1:12,labels=substr(month.names,1,1))
for (i in 1:15) lines(1:12,max.data[i,],col=i)
dev.off()

