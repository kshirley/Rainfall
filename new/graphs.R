#Comes from the graphs in R_code_tobit_real3.R and the graphs in comparison_plots2.R
#After running the preferred baseline file, this file should work to create all graphs using the resulting posterior samples.

#There are several sets of graphs: 
#1. Trace plots of parameters
#2. Onset metric
#3. Consecutive dry days metric
#4. Probability/Mean/Sd/Max Posterior Plots
#5. Historical Plots
#6. Comparison Plots: Posterior to Historical with line plots
#7. Comparison Plots: Posterior to Historical with error bars
#

# Clear the workspace:
rm(list=ls())
setwd("~/Git/Rainfall/")

path<-"/Users/kathrynvasilaky/Documents/OneDrive/IRI/RainfallSimulation/Rainfall/Rainfall/"
#source(paste(path,"/Rcode_tobit_mcmc_functions.R",sep=""))
#Pulls all the functions necessary for plotting
source("Rcode_tobit_mcmc_functions.R")
setwd("/Users/kathrynvasilaky/Documents/OneDrive/IRI/RainfallSimulation/Rainfall/Rainfall")



###############################
#Load last MCMC
###############################

# load the input data and assign them to the global namespace:
load("input_data.RData")
for (i in 1:length(input.data)) assign(names(input.data)[i], input.data[[i]])


#Regular Posterior
#gibbs_out_20150326_G5k.RData
#gibbs_out_20150616_G20k.RData
#gibbs_out_20150326_G5k_HS1.RData
#gibbs_out_20150514_G10k.RData
#Last 3 years of Hager Salem as NA:
load(file="gibbs_out_20150616_G20k.RData")
for (i in 1:length(gibbs.list)) assign(names(gibbs.list)[i],gibbs.list[[i]])


############################
#Graph Posterior Parameters
############################
#R vs W
#Problem is that W is not stored
pdf(file=paste(path,"fig_Y_vs_W.pdf",sep=""),width=7,height=7)
par(mfrow=c(1,1))
for (s in 1:S){
  for (j in 1:J[s]){
    plot(Y[[s]][j,Y[[s]][j,]>0])
    points(W[[s]][j,is.na(Y[[s]][j,])],col=2)
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
tp(lambda.gibbs,burn=1500)
#abline(h=true.values$lambda,col=4)
dev.off()

# alpha (t distribution degrees of freedom)
pdf(file=paste(path,"fig_alpha.pdf",sep=""),width=7,height=7)
par(mfrow=c(1,1))
tp(alpha.gibbs)
tp(alpha.gibbs,burn=1500)
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
for (i in 1:12) se.nino[i] <- sd(as.numeric(mu.gibbs[,1501:G,i+11]))
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


####Metrics: Onset, Consecutive Rainfall, Prob/Mean/Max/SD Rainfall#####

##############
#Onset Metric
#############
#Measured as the first day in July on which cumulative July rainfall exceeds 20 m
#Self Contained vs the rest of the metrics which compute the data above and the posterior data below

#20 mm or more over 3 consecutive rainy days after July 1st
# Measure onset in the observed data:
# where onset date = the day in April on which cumulative April rainfall exceeded 20 mm
july.obs <- as.list(rep(NA, S))
july.sum <- as.list(rep(NA, S))
for (s in 1:S) {
  july.obs[[s]] <- matrix(NA, 19, J[s])
  july.sum[[s]] <- matrix(NA, 19, J[s])
  for (l in 1:J[s]) {
    for (year in 1:19) {
      sel <- as.numeric(substr(date.string, 1, 4)) == 1991 + year & months(date.string) == "July"
      july.obs[[s]][year, l] <- sum(!na.mat[[s]][l, sel])
      july.sum[[s]][year, l] <- sum(Y[[s]][l, sel])
    }    
  }
}

#####Real data#####
# for July with 30 observed days of rainfall, measure day on which cumulative July sum surpassed 20mm:
onset.day <- as.list(rep(NA, S))
#days 31+31-2-2
tmp<-c(rep(NA,30))
for (s in 1:S) {
  onset.day[[s]] <- matrix(NA, 19, J[s])
  for (l in 1:J[s]) {
    for (year in 1:19) {
      sel <- as.numeric(substr(date.string, 1, 4)) == 1991 + year & months(date.string) == "July" 
      for (j in 3:30) {
        #tmp <- cumsum(R[[s]][l, sel])
        #this is the statement to say that the last two days have to have been positive
        if ( Y[[s]][l,sel][j-1]>0 & Y[[s]][l,sel][j-2]>0 & !is.na(Y[[s]][l,sel][j-1]) & !is.na(Y[[s]][l,sel][j-2]) ) {
          tmp[i]=Y[[s]][l,sel][j]+Y[[s]][l,sel][j-1]+Y[[s]][l,sel][j-2] 
        }
        else {tmp[j]=NA}
        onset.day[[s]][year, l] <- ifelse(!is.na(min(which(tmp > 20))), min(which(tmp > 20))+2, NA)
      }  
    }
  }
}


#####Simulated Data for Onset#####
# Number of posterior predictive simulations
N.sim <- 150

# Set up place to store posterior predictive simulations:
pp.onset.date <- as.list(rep(NA, S))
for (s in 1:S) {
  pp.onset.date[[s]] <- array(NA, dim=c(19, J[s], N.sim))
}


# Spread out (thin) across chains and iterations:
k.vec <- rep(1:3, 50)
g.vec <- rep(round(seq(1001, 2000, length=50)), each=3)

# Set up year selection and store to access within loop:
sel.year <- as.list(rep(NA, 19))
for (i in 1:19) {
  sel.year[[i]] <- which(as.numeric(substr(date.string, 1, 4)) == 1991 + i & months(date.string) == "July")  
}

tmp<-c(rep(NA,30))
# Do the posterior predictive simulations:
#KV-With Thinning
t1 <- Sys.time()
for (i in 1:N.sim) {
  if (i %% 10 == 0) print(i)
  # Collect Sigma into a list:
  Sigma.sim <- as.list(rep(NA, S))
  for (s in 1:S) Sigma.sim[[s]] <- Sigma.gibbs[[s]][k.vec[i], g.vec[i], , ]
  
  # Simulate rainfall:
  W.new <- sim.W(alpha=alpha.gibbs[k.vec[i], g.vec[i]], beta=beta.gibbs[k.vec[i], g.vec[i], , ], 
                 lambda=lambda.gibbs[k.vec[i], g.vec[i]], 
                 tau=tau.gibbs[k.vec[i], g.vec[i]], beta.arc=beta.arc.gibbs[k.vec[i], g.vec[i], ],
                 Sigma=Sigma.sim, X=X, X.arc=X.arc, na.preserve=TRUE, na.mat=na.mat)
  for (s in 1:S) {
    for (l in 1:J[s]) {
      for (year in 1:19) {
        #from the 3rd day to the 30th
        for (j in 3:30) {
          #tmp <- cumsum(R[[s]][l, sel])
          #this is the statement to say that the last two days have to have been positive
          if ( W.new[[s]][l,sel.year[[year]]][j-1]>0 & W.new[[s]][l,sel.year[[year]]][j-2]>0 & !is.na(W.new[[s]][l,sel.year[[year]]][j-1]) & !is.na(W.new[[s]][l,sel.year[[year]]][j-2]) ) {
            #fill the temp matrix of 30 days with consecutive rainfall if the previous 2 days had rainfall  
            tmp[j]=W.new[[s]][l,sel.year[[year]]][j]+W.new[[s]][l,sel.year[[year]]][j-1]+W.new[[s]][l,sel.year[[year]]][j-2] 
          }
          #else {tmp[j]=NA}
          #pull out the day +2 on wihch rainfall was greater than 20
          pp.onset.date[[s]][year, l, i] <- ifelse(!is.na(min(which(tmp > 20))), min(which(tmp > 20)+2), NA)    
        }
      }
    }
  }
}

t2 <- Sys.time()
t2 - t1
# 11 seconds for N.sim = 150 ! (yay.)

# Let's make ~ 170 posterior predictive plots!
pdf(file="Onset_postpred.pdf", width=6, height=6)
par(mfrow=c(2, 2))
for (s in 1:S) {
  for (year in 1:19) {
    for (l in 1:J[s]) {
      if (july.obs[[s]][year, l] == 30) {
        #pdf(paste("Onset Rainfall", i,".pdf",sep=""))
        hist(pp.onset.date[[s]][year, l, ], breaks=0:30, main="", las=1, xlab="Day in July")
        legend("topleft", legend=paste(sum(is.na(pp.onset.date[[s]][year, l, ])), "/150 missing", sep=""), inset=0.02)
        title(main=paste(rownames(Y[[s]])[l], year + 1991))
        abline(v = onset.day[[s]][year, l], col="orange", lwd=2)
        abline(v = mean(pp.onset.date[[s]][year, l, ], na.rm=TRUE), col=2, lty=2)
      }
    }
  }
}
dev.off()





#################################
#Consecutive Dry Days Collection
#################################
J.sum <- sum(J)
# First code block inserted
consec_zeros_count <- function(temp) {
  #This function takes in a list of rainfall counts and computes
  #the length of each set of contiguous zeros.
  temp2 <- temp
  temp2[temp!=0] <- 0
  temp2[is.na(temp)] <- 0
  temp2[temp==0] <- 1
  rle_out <- rle(temp2)
  len <- rle_out$lengths
  v <- rle_out$values
  consec_zeros <- len[v==1]
}

#bins for histograms
edges <- seq(0,100)
#add a bin for everything larger than 100 days
edges <- append(edges,10000)
#kenny's bayesian thinning, select index in every 50 after burn in
g.vec <- round(seq(3001, 5000, length=50))
#creat empty data structure filled with na's for consec hist
zero_count_data.gibbs <- array(NA,dim = c(K, length(g.vec), J.sum, length(edges)-1))
#end insertion

#Generate W.new data, from posterior distributions
#KV-K is number of chains, G is number of runs per chain, J.sum is the total number of sereis across locations
#KV-right now it generates 5000 runs
cumsum.data.gibbs<- array(NA, dim=c(K,length(g.vec),J.sum,19)) #annual total rainfall

for (k in 1:K){
  print( k)
  for (g in 10001:20000){
    if (g%%20==0) { print(g) }
    
    # Collect Sigma into a list:
    Sigma.sim <- as.list(rep(NA,S))
    for (s in 1:S) Sigma.sim[[s]] <- Sigma.gibbs[[s]][k,g,,]
    # Simulate rainfall:
    W.new <- sim.W(alpha=alpha.gibbs[k,g],beta=beta.gibbs[k,g,,],lambda=lambda.gibbs[k,g],tau=tau.gibbs[k,g],beta.arc=beta.arc.gibbs[k,g,],
                   Sigma=Sigma.sim,X=X,X.arc=X.arc,na.preserve=TRUE,na.mat=na.mat)
    
    #This is the second of 2 code blocks that inserted
    #check to make sure the current g is at least as big as the smallest g.vec
    if (g >= min(g.vec)) {
      #only capture the correct slices
      #if g is an element in the gvec, or thinned series
      if (g %in% g.vec) {
        #need to find where we are at in g.vec
        # valid answers should be from 1 to length(g.vec)
        local_g <- which(g.vec==g)
        for (s in 1:S){
          for (j in 1:J[s]){
            #create vector of dry spells
            dry_days <- consec_zeros_count(W.new[[s]][j,])
            #count frequency of dryspells stored in sim data, where bins are constant across trials, but we only pull out 50 trials
            #think a matrix 15 by 6779, repeated 5000 times, but we only pull out 50 of those trials. 
            #compress those 6779 days into 100 elemenet vector
            dry_days_hist <- hist(dry_days, breaks=edges)
            zero_count_data.gibbs[k,local_g, c(0,cumsum(J))[s]+j,] <- dry_days_hist$counts
            
            #cumsum.data.gibbs[3,50,10,19] <- round(sum(W.new[[s]][j,year.vec==y],na.rm=TRUE))
            
            #end insertion
            #for (m in 1:12){
            for (y in year.vec){
              #new for interannual variability
              cumsum.data.gibbs[k,local_g,c(0,cumsum(J))[s]+j,y] <- round(sum(W.new[[s]][j,year.vec==y],na.rm=TRUE))
              print(round(sum(W.new[[s]][j,year.vec==y],na.rm=TRUE)))
              print( k)
            } 
          }
        }
      }
    }
  }
}



##########################
#Posterior Data Collection
##########################
#Probability Wet Days
#Mean Wet Days
#Variance Wet Days
#Max Wet Days
####################
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
#save(post.list,file=paste(path,"post_list_G5000_092614.RData",sep=""))
#load(file=paste(path,"post_list_G5000_092614.RData",sep=""))
#for (i in 1:length(post.list)) assign(names(post.list)[i],post.list[[i]])


#save all the data sets for future use, so you don't have to rerun for 5 hours
#save.image("~/SkyDrive/IRI/RainfallSimulation/R/Rainfall/post_sim_output092614.RData")



#############################
#Historical Data Collection
#############################
# compute the mean monthly rainfall:
#month.vec <- match(months(date.string),month.names)
#year.vec <- match(year.names,year.names)
J.sum <- sum(J)


#KV-Now collect the historical???? data to plot
# From the data:
p.wet.data <- matrix(NA,J.sum,12)
mean.data <- matrix(NA,J.sum,12)
sd.data <- matrix(NA,J.sum,12)
max.data <- matrix(NA,J.sum,12)
cumsum.data<-matrix(NA, J.sum, 19)

#wet.vec<-matrix(NA, J.sum, 19)

for (s in 1:S){
  for (j in 1:J[s]){
    for (m in 1:12){
      for (y in year.vec){
        cumsum.data[c(0,cumsum(J))[s]+j,y] <- round(sum(Y[[s]][j,year.vec==y],na.rm=TRUE))
    
        #print(c(0,cumsum(J))[s]+j)
        #R is the list of 
        p.wet.data[c(0,cumsum(J))[s]+j,m] <- sum(Y[[s]][j,month.vec==m]>0,na.rm=TRUE)/sum(!is.na(Y[[s]][j,month.vec==m]))
        wet.vec <- Y[[s]][j,month.vec==m & Y[[s]][j,]>0]    
        mean.data[c(0,cumsum(J))[s]+j,m] <- mean(wet.vec,na.rm=TRUE)
        sd.data[c(0,cumsum(J))[s]+j,m] <- sd(wet.vec,na.rm=TRUE)
        max.data[c(0,cumsum(J))[s]+j,m] <- max(wet.vec,na.rm=TRUE)
        #print(paste(s,j,m,length(wet.vec),sum(!is.na(wet.vec))),sep=",")
      }
    }
  }
}

#this pulls out the consec dry days across 6779 days per site and histograms into 100 bins
consec.data<- matrix(data=NA, nrow=J.sum, ncol=length(edges)-1)
for (s in 1:S){
  for (j in 1:J[s]){
    wet.vec2 <- Y[[s]][j,]
    dry_days<-consec_zeros_count(wet.vec2)
    dry_days_hist <- hist(dry_days, breaks=edges)
    #this cumsum(J))[s]+j, is a counter going from 1-15th row
    consec.data[c(0,cumsum(J))[s]+j,] <- dry_days_hist$counts
  }
}

dev.off()

# Replace the Infs with NAs else you'll get an error for plots
max.data[max.data==-Inf] <- NA

#KV-Plot the actual original data:
par(mfrow=c(1,1))


# probability of rain:

pdf(file=paste(path,"PercWetDays_Historical.pdf",sep=""),height=6,width=9)
plot(-10,-10,xlim=c(1,12),ylim=c(0,1),las=1,xlab="Month",ylab="% Wet Days",xaxt="n",main="% of wet days, HISTORICAL")
axis(1,at=1:12,labels=substr(month.names,1,1))
for (i in 1:15) lines(1:12,p.wet.data[i,],col=i)
dev.off()


# mean wet day rainfall:
pdf(file=paste(path,"MeanWetDays_Historical.pdf",sep=""),height=6,width=9)
plot(-10,-10,xlim=c(1,12),ylim=c(0,max(mean.data,na.rm=TRUE)),las=1,xlab="Month",ylab="Mean (mm)",xaxt="n",main="Mean wet day rainfall,HISTORICAL")
axis(1,at=1:12,labels=substr(month.names,1,1))
for (i in 1:15) lines(1:12,mean.data[i,],col=i)
dev.off()


# sd wet day rainfall:
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
dev.off()

###########################
#Posterior Histogram Plots
###########################
# Plot these posterior predictive checks:
pb <- seq(1001,2000,by=5)
pdf(file=paste(path,"fig_tobit_sd.gibbs.pdf",sep=""),width=10,height=8)
par(mfrow=c(3,4))
for (i in 1:J.sum){
  for (m in 1:12){
    #g.out <- p.wet.gibbs[,pb,i,m]; d.out <- p.wet.data[i,m]
    #lab <- "P(Wet)"
    #g.out <- mean.gibbs[,pb,i,m]; d.out <- mean.data[i,m]
    #lab <- "Mean wed days (mm)"
    g.out <- sd.gibbs[,pb,i,m]; d.out <- sd.data[i,m]
    lab <- "Sd wet days (mm)"
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


# Plot some values of W against the observed values (R):
pdf(file=paste(path,"fig_ppred_WvsObservedR.pdf",sep=""),width=10,height=8)
par(mfrow=c(2,1))
plot(W.new[[1]][1,na.mat[[1]][1,]==0])
plot(Y[[1]][1,na.mat[[1]][1,]==0],col=2)
dev.off()


#Trace plot of covariance matrix elements:
pdf(file=paste(path,"fig_ppred_cov_matrix.pdf",sep=""),width=10,height=8)
s <- 1
par(mfrow=c(2,2))
tp(Sigma.gibbs[[s]][,,1,1],las=1)
tp(Sigma.gibbs[[s]][,,1,2],las=1)
tp(Sigma.gibbs[[s]][,,2,1],las=1)
tp(Sigma.gibbs[[s]][,,2,2],las=1)
dev.off()



#KV-More posterior predictive checks:

####################
#Posterior Pvalues
####################
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

#KV-Plot the distribution of the pvalues for all the moments from the posterior-simulated W's
pdf(file=paste(path,"fig_pvals_hist.pdf",sep=""),width=10,height=8)
stat.names <- c("P(wet)","Mean","Sd","Max")
par(mfrow=c(2,2))
for (i in 1:4) hist(p.val[i,,],xlim=c(0,1),breaks=seq(0,1,0.1),main=stat.names[i])
dev.off()

############################
#Posterior Plots by Series
############################

#KV plot comparable posterior graphs to historical
post <- p.wet.gibbs[1,10001:20000,1:15,1:12]
post_p <- colMeans(post)

pdf(file=paste(path,"PercWetDays_Posterior.pdf",sep=""),height=6,width=9)
plot(-10,-10,xlim=c(1,12),ylim=c(0,1),las=1,xlab="Month",ylab="% Wet Days",xaxt="n",main="% of wet days, Posterior")
axis(1,at=1:12,labels=substr(month.names,1,1))
for (i in 1:15) lines(1:12,post_p[i,],col=i)
dev.off()

# mean wet day rainfall:
post <- mean.gibbs[1,10001:20000,1:15,1:12]
post_m <- colMeans(post)

pdf(file=paste(path,"MeanWetDays_Posterior.pdf",sep=""),height=6,width=9)
plot(-10,-10,xlim=c(1,12),ylim=c(0,max(post_m,na.rm=TRUE)),las=1,xlab="Month",ylab="Mean (mm)",xaxt="n",main="Mean wet day rainfall,Posterior")
axis(1,at=1:12,labels=substr(month.names,1,1))
for (i in 1:15) lines(1:12,post_m[i,],col=i)
dev.off()


# sd wet day rainfall:
post <- sd.gibbs[1,10001:20000,1:15,1:12]
post_s <- colMeans(post)

pdf(file=paste(path,"SdWetDays_Posterior.pdf",sep=""),height=6,width=9)
plot(-10,-10,xlim=c(1,12),ylim=c(0,max(post_s,na.rm=TRUE)),las=1,xlab="Month",ylab="Sd (mm)",xaxt="n",main="Sd of wet day rainfall,Posterior")
axis(1,at=1:12,labels=substr(month.names,1,1))
for (i in 1:15) lines(1:12,post_s[i,],col=i)
dev.off()


# max wet day rainfall:
post <- max.gibbs[1,10001:20000,1:15,1:12]
post_ma <- colMeans(post)

pdf(file=paste(path,"MaxWetDays_Posterior.pdf",sep=""),height=6,width=9)
plot(-10,-10,xlim=c(1,12),ylim=c(0,max(post_ma,na.rm=TRUE)),las=1,xlab="Month",ylab="Max (mm)",xaxt="n",main="Max of wet day rainfall,Posterior")
axis(1,at=1:12,labels=substr(month.names,1,1))
for (i in 1:15) lines(1:12,post_ma[i,],col=i)
dev.off()


##############################
#Posterior Plots as Error Bars
#from comparision plots1.R
##############################


#KV-this squishes the data to 15 by 12, by averaging over the last 1000 MC runs for the first chain
#for each location (15) and for eahc month (12)
#p.wet.gibbs[K chain,G observations,J series,m months]
#p.wet.data is from the historical original data

#comparision prob plot for each location of posterior bar plot with the historical data
pdf(file=paste(path,"Posterior_pwet_series.pdf",sep=""),width=6,height=9)
par(mfrow=c(3,3))
for (i in 1:J.sum){
  post <- p.wet.gibbs[1,1:50,i,1:12]
  #pdf(paste("Posterior_pwet_series", i,".pdf",sep=""))
  boxplot(post, xaxt = "n", xlim=c(1,12), ylim=c(0,1), main=i)
  axis(1,at=1:12,labels=substr(month.names,1,1), xlim=c(1,12), ylim=c(0,1))
  lines(1:12, p.wet.data[i,], col="dark red")
  #dev.off()
}
dev.off()

#comparision mean plot for each location of posterior bar plot with the historical data
pdf(file=paste(path,"Posterior_mean_series.pdf",sep=""),width=6,height=9)
par(mfrow=c(3,3))
for (i in 1:J.sum){
  post <- mean.gibbs[1,1:50,i,1:12]
  #pdf(paste("Posterior_mean_series", i,".pdf",sep=""))
  boxplot(post, xaxt = "n", xlim=c(1,12), ylim=c(0,30),main=i)
  axis(1,at=1:12,labels=substr(month.names,1,1), xlim=c(1,12), ylim=c(0,30))
  lines(1:12, mean.data[i,],col="dark red")
  #dev.off()
}
dev.off()


#comparision mean plot for each location of posterior bar plot with the historical data
pdf(file=paste(path,"Posterior_max_series.pdf",sep=""),width=6,height=9)
par(mfrow=c(3,3))
for (i in 1:J.sum){
  post <- max.gibbs[1,1:50,i,1:12]
  #pdf(paste("Posterior_max_series", i,".pdf",sep=""))
  boxplot(post, xaxt = "n", xlim=c(1,12), ylim=c(0,150),main=i)
  axis(1,at=1:12,labels=substr(month.names,1,1), xlim=c(1,12),ylim=c(0,150))
  lines(1:12, max.data[i,],col="dark red")
  #dev.off()
}
dev.off()

#comparison of consecutive dry days for simulated and historical data
pdf(file=paste(path,"Consec_Data_Series.pdf",sep=""),width=6,height=9)
par(mfrow=c(3,3))
for (i in 1:J.sum){
  #for log scale use log10="y", just 1 chain here
  #pdf(paste("Consec_Data_Series", i,".pdf",sep=""))
  boxplot(zero_count_data.gibbs[1,1:50,i,], xlim=range(c(0,30)),ylim=range(c(0,250)),main= i)
  par(new=TRUE)
  #axes=FALSE drops one axis, type="s" is step plot
  plot(consec.data[i,], xlim=range(c(0,30)), ylim=range(c(0,250)), ylab="Frequency in 6779 days", xlab="count dry days", axes=FALSE,col = "dark red", type="s")
  #dev.off()
}
dev.off()


#comparison of total annual rainfall for simulated and historical data
pdf(file=paste(path,"Total_Annual_Rainfall.pdf",sep=""),width=6,height=9)
par(mfrow=c(2,2))
for (i in 1:J.sum-1){
  #for log scale use log10="y", just 1 chain here
  #pdf(paste("Total_Annual_Rainfall", i,".pdf",sep=""))
  boxplot(cumsum.data.gibbs[1,1:5,i,], xlim=range(c(0,19)),ylim=range(c(250,1050)), axes=TRUE,main= i)
  par(new=TRUE)
  #axes=FALSE drops one axis, type="s" is step plot
  #axis(1,at=1:19,labels=substr(year.names,1,1), xlim=c(0,19),ylim=c(250,550))
  lines(cumsum.data[i,], col="dark red")
  #xlim=range(c(0,19)), ylim=range(250,550), ylab="Total Rainfall", xlab="Years", axes=FALSE,col = "dark red", type="s")
  #dev.off()
}
dev.off()


#comparison of interannual rainfall for simulated and historical data
#cumsum.data.gibbs[k=runs,g=splice,sj=series,y=year]
#fix series to first run, gen sd's for each run across years
sd.data.gibbs<- array(NA, dim=c(length(g.vec),J.sum)) #sd total rainfall across years by series

for (g in 1:50){
  for (j in 1:J.sum){
    print(sd(cumsum.data.gibbs[,g,j,], na.rm=TRUE))
    sd.data.gibbs[g,j]<-sd(cumsum.data.gibbs[,g,j,], na.rm=TRUE)   
  }}

sd.data<-array(NA, 15)
for (i in 1:J.sum) {
  sd.data[i]<-sd(cumsum.data[i,])
}

pdf(paste("Interannual",".pdf",sep=""))
boxplot(sd.data.gibbs[,],xlab="Site and Series", ylab="Simulated and Historical (RED) Interannual Variability by Site")
par(new=TRUE)
lines(sd.data,col="dark red")
dev.off()

