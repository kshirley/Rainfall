
#This file sets the starting points for the within site covariance Sigma.start
#It uses the ending values form a previous run to set the mean starting points for chain 1, 4 sds above and below for chains 2 and 3

setwd("~kathrynvasilaky/Documents/OneDrive/IRI/RainfallSimulation/Rainfall/Rainfall")
path<-"kathrynvasilaky/Documents/OneDrive/IRI/RainfallSimulation/Rainfall/Rainfall"


S <- 6
K <- 3
Sigma.start <- vector(mode = "list", length = S)
for (s in 1:S) {
  Sigma.start[[s]] <- array(0, dim=c(K, J[s], J[s]))
}

load(file = "gibbs_out_20150514_G10k.RData")
for (i in 1:length(gibbs.list)) assign(names(gibbs.list)[i], gibbs.list[[i]])

#Loop across first 5 sites
for (s in 1:5) {
# site, first chain
  #Sigma.start[[site]][chain, row , column ]
Sigma.start[[s]][1, 1 ,2 ] <- mean(Sigma.gibbs[[s]][1, , 1, 1])
Sigma.start[[s]][1, 2 ,1 ] <- mean(Sigma.gibbs[[s]][1, , 2, 1])
Sigma.start[[s]][1, 1 ,1 ] <- mean(Sigma.gibbs[[s]][1, , 1, 1])
Sigma.start[[s]][1, 2 ,2 ] <- mean(Sigma.gibbs[[s]][1, , 2, 2])



# site, second chain, set it 4 sds above the mean
Sigma.start[[s]][2, 1 ,2 ] <- mean(Sigma.gibbs[[s]][2, , 1, 1]) + 4*sd(Sigma.gibbs[[s]][2, , 2, 1])
Sigma.start[[s]][2, 2 ,1 ] <- mean(Sigma.gibbs[[s]][2, , 2, 1]) + 4*sd(Sigma.gibbs[[s]][2, , 2, 1])
Sigma.start[[s]][2, 1 ,1 ] <- mean(Sigma.gibbs[[s]][2, , 1, 1]) + 4*sd(Sigma.gibbs[[s]][2, , 1, 1])
Sigma.start[[s]][2, 2 ,2 ] <- mean(Sigma.gibbs[[s]][2, , 2, 2]) + 4*sd(Sigma.gibbs[[s]][2, , 2, 2])

# site, third chain, set it 4 sds below the mean
Sigma.start[[s]][3, 1 ,2 ] <- mean(Sigma.gibbs[[s]][3, , 2, 1]) - 4*sd(Sigma.gibbs[[s]][3, , 2, 1])
Sigma.start[[s]][3, 2 ,1 ] <- mean(Sigma.gibbs[[s]][3, , 1, 1]) - 4*sd(Sigma.gibbs[[s]][3, , 2, 1])
Sigma.start[[s]][3, 1 ,1 ] <- mean(Sigma.gibbs[[s]][3, , 1, 1]) - 4*sd(Sigma.gibbs[[s]][3, , 1, 1])
Sigma.start[[s]][3, 2 ,2 ] <- mean(Sigma.gibbs[[s]][3, , 2, 2]) - 4*sd(Sigma.gibbs[[s]][3, , 2, 2])

}


#Loop across rows and columns for last site
#Sigma.start[[site]][chain, row , column ]



for (i in 1:5) {
  for (j in 1:5) {
    #First chain, 6th site
    Sigma.start[[6]][1, i ,j ] <- mean(Sigma.gibbs[[6]][1, , i, j])    
  }
}

for (i in 1:5) {
  for (j in 1:5) {
    #Second chain, 6th site
    Sigma.start[[6]][2, i ,j ] <- mean(Sigma.gibbs[[6]][2, , i, j]) + 4*sd(Sigma.gibbs[[6]][2, , i, j]) 
  }
}



#this site has sparse data and does not like 4 sds below the mean so it's at 1 sd below. 
for (i in 1:5) {
  for (j in 1:5) {
    #Third chain, 6th site
    Sigma.start[[6]][3, i ,j ] <- mean(Sigma.gibbs[[6]][3, , i, j]) - 1*sd(Sigma.gibbs[[6]][3, , i, j])   
  }
}



