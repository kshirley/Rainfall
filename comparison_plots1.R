#New Plots for paper

#load("~/SkyDrive/IRI/RainfallSimulation/R/Rainfall/post_sim_output042914.RData")

load("~/SkyDrive/IRI/RainfallSimulation/Rainfall/gibbs_out_NA11172014_G5000.RData")

'''
K-3 chains
S-6 sites
G-5000 MC observations
J-15 series, 2 for 4,and 5 for last 
J.sum-15
M-12 months
'''


#IN R_code_tobit_real2.R :
  
#  line 560 is the thinning of the series for a different metric

#line 581 begins the consec dry day metric 

#What's need is a plot of the posterior distribution of the consec dry days (with error bars) and the historical plot of consec dry days drawn through. 



#KV-this squishes the data to 15 by 12, by averaging over the last 1000 MC runs for the first chain
#for each location (15) and for eahc month (12)
#p.wet.gibbs[K chain,G observations,J series,m months]
#p.wet.data is from the historical original data

post <- p.wet.gibbs[1,4001:5000,1:15,1:12]
post_p <- colMeans(post)

#comparision prob plot for each location of posterior bar plot with the historical data
for (i in 1:J.sum){
  post <- p.wet.gibbs[3,4001:5000,i,1:12]
  pdf(paste("Posterior_pwet_series", i,".pdf",sep=""))
  boxplot(post, xaxt = "n", xlim=c(1,12), ylim=c(0,1))
  axis(1,at=1:12,labels=substr(month.names,1,1), xlim=c(1,12), ylim=c(0,1))
  lines(1:12, p.wet.data[i,])
  dev.off()
}


#comparision mean plot for each location of posterior bar plot with the historical data
for (i in 1:J.sum){
  post <- mean.gibbs[3,4001:5000,i,1:12]
  pdf(paste("Posterior_mean_series", i,".pdf",sep=""))
  boxplot(post, xaxt = "n", xlim=c(1,12), ylim=c(0,30))
  axis(1,at=1:12,labels=substr(month.names,1,1), xlim=c(1,12), ylim=c(0,30))
  lines(1:12, mean.data[i,])
  dev.off()
}

#comparision mean plot for each location of posterior bar plot with the historical data
for (i in 1:J.sum){
  post <- max.gibbs[3,4001:5000,i,1:12]
  pdf(paste("Posterior_max_series", i,".pdf",sep=""))
  boxplot(post, xaxt = "n", xlim=c(1,12), ylim=c(0,150))
  axis(1,at=1:12,labels=substr(month.names,1,1), xlim=c(1,12),ylim=c(0,150))
  lines(1:12, max.data[i,])
  dev.off()
}

#comparison of consecutive dry days for simulated and historical data
for (i in 1:J.sum){
  #for log scale use log10="y", just 1 chain here
  pdf(paste("Consec_Data_Series", i,".pdf",sep=""))
  boxplot(zero_count_data.gibbs[1,1:50,i,], xlim=range(c(0,30)),ylim=range(c(0,250)))
  par(new=TRUE)
  #axes=FALSE drops one axis, type="s" is step plot
  plot(consec.data[i,], xlim=range(c(0,30)), ylim=range(c(0,250)), ylab="Frequency in 6779 days", xlab="count dry days", axes=FALSE,col = "dark red", type="s")
  dev.off()
}


#comparison of total annual rainfall for simulated and historical data
for (i in 1:J.sum-1){
  #for log scale use log10="y", just 1 chain here
  pdf(paste("Total Annual Rainfall", i,".pdf",sep=""))
  boxplot(cumsum.data.gibbs[1,1:5,i,], xlim=range(c(0,19)),ylim=range(c(250,1050)), axes=TRUE)
  par(new=TRUE)
  #axes=FALSE drops one axis, type="s" is step plot
  #axis(1,at=1:19,labels=substr(year.names,1,1), xlim=c(0,19),ylim=c(250,550))
  lines(cumsum.data[i,], col="dark red")
        #xlim=range(c(0,19)), ylim=range(250,550), ylab="Total Rainfall", xlab="Years", axes=FALSE,col = "dark red", type="s")
    dev.off()
}

#comparison of interannual rainfall for simulated and historical data
#cumsum.data.gibbs[k=runs,g=splice,sj=series,y=year]
#fix series to first run, gen sd's for each run across years
sd.data.gibbs<- array(NA, dim=c(length(g.vec),J.sum)) #sd total rainfall across years by series

for (g in 1:5){
    for (j in 1:J.sum){
    print(sd(cumsum.data.gibbs[,g,j,], na.rm=TRUE))
    sd.data.gibbs[g,j]<-sd(cumsum.data.gibbs[,g,j,], na.rm=TRUE)
    
}}

sd.data<-array(NA, 15)
for (i in 1:J.sum) {
  sd.data[i]<-sd(cumsum.data[i,])
}

boxplot(sd.data.gibbs[,],xlab="Site and Series", ylab="Simulated and Historical (RED) Interannual Variability by Site")
par(new=TRUE)
lines(sd.data,col="dark red")
pdf(paste("Interannual",".pdf",sep=""))
dev.off()
      






#comparision plot for each location of posterior bar plot with the historical data
#I want to sum data to get cumulative rainfall 
#What's complicated is that I want to count the frequency of consec dry days for atleast 100 different runs, 
#that way I have 100 data points for each day in a series, so I can get bar charts
#But....those vectors of frequencies will be of different lengths in each run

#http://stackoverflow.com/questions/5012516/count-how-many-consecutive-values-are-true
#http://stackoverflow.com/questions/1502910/how-can-i-count-runs-in-a-sequence


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

freq<- matrix(NA,S,J.sum)
freq=list()

for (s in 1:S){
  freq[[s]] <- array(NA,dim=c(J[s], 1000))
  for (j in 1:J[s]){
      temp<-consec_zeros_count(W.new[[s]][j,])
      print(length(temp))
      freq[[s]][j, 1:length(temp)] <- temp
}}


#This double loop is bad, ignore!
for (s in 1:S){
  for (j in 1:J[s]){
    data <- as.vector(zero_count_data.gibbs[1:3,1:100,c(0,cumsum(J))[s]+j,])
    temp <- consec_zeros_count(data)
}}


temp <- array(NA,dim=c(100,1000))
for (s in 1:S){
  for (j in 1:J[s]){
    
    for (t in 1:100){
      data <- zero_count_data.gibbs[k,t,c(0,cumsum(J))[s]+j,]
      temp[consec_zeros_count(data)
}}}



# Spread out (thin) across chains and iterations:
k.vec <- rep(1:3, 50)
g.vec <- round(seq(3001, 5000, length=50))

edges <- seq(0,100)
edges <- append(edges,10000)
bob <- hist(temp, breaks=edges)


for (i in seq(1000,5000)){
  if (i %in% g.vec) {
    ndx<-which(g.vec==i)
    print(i)
    print(ndx)
}}


