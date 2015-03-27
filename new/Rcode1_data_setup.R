# clear workspace:
rm(list=ls())

# set working directory:
setwd("~/Git/Rainfall/")

#setwd("~/SkyDrive/IRI/RainfallSimulation/Rainfall")
#path<-"~/SkyDrive/IRI/RainfallSimulation/Rainfall/"



### 1. The observed rainfall data ### 

# Load in 15 time series for Ethiopia:
load("ethiopia_full_data.RData")

# this is a 15 x 18094 data.frame, where the rows correspond
# to rainfall time series (they are labeled) and the columns
# are days starting with 1/1/1961, ending 7/28/2010, 
# excluding Feb 29 from all leap years: (1964, 1968, 1972, 1976, 
# 1980, 1984, 1988, 1992, 1996, 2000, 2004, 2008)
# This means 12 days are excluded.

# check length of series:
length(seq(as.Date("1961-01-01"), as.Date("2010-07-28"), by = 1))
# 18106 days, including leap year days

# our data:
dim(data)[2]
# 18094, 12 fewer days. Good, this matches the number of leap years

# Re order the rows, to group by location:
data <- data[c(6, 1, 7, 2, 8, 3, 9, 4, 10, 5, 11, 12, 13, 14, 15), ]

# set the data to 1992 onward, so get rid of the first 31 years:
data <- data[, (365*31 + 1):dim(data)[2]] # 1992 onward, all 15 sites

# get location names:
site.names <- unlist(strsplit(rownames(data), "_")[seq(1, by=2, length=6)])[seq(1, by=2, length=6)]

# Define number of days in the time series:
N <- dim(data)[2]

# number of locations:
S <- 6

# number of series for each location
J <- c(2, 2, 2, 2, 2, 5)

# Create the datestring beginning at Jan 1 1992:
# Manually add 5 leap days:
date.string <- as.Date(1:(N + 5) - 1, origin = "1992-01-01")

# Now define the indices of those 5 leap days to delete:
leap.year.index <- 60 + c(0, 365*4 + 1, 365*8 + 2, 365*12 + 3, 365*16 + 4)
date.string <- date.string[-leap.year.index]

# Get month names and year names, and create vectors to store their value for each day:
month.names <- unique(months(date.string))
year.names <- unique(format(date.string, "%Y"))
month.vec <- match(months(date.string), month.names)
year.vec <- match(year.names, year.names)


# KV: IF YOU'RE CUTTING OFF A SERIES DO IT HERE
# Make all series for last three years of Hagersalam, NA, for cross validation exericse
# This is 5476 to 6779 is 2007 to 2010
#data[1, 5476:6779]<-NA

# Put the observed rainfall in a list with one element per location
# where each list element is a matrix:
Y <- as.list(rep(NA, S))
for (s in 1:S) Y[[s]] <- data[1:J[s] + c(0, cumsum(J))[s], ]




### 2. The lat-long location data ### 
ll <- read.table("lat_long_details.csv", sep=",", header=TRUE)
ll <- ll[c(3, 5, 4, 1, 2, 6), ]  # reorder to match site.names ordering
d.mat <- matrix(NA, S, S)
for (i in 1:S){
  for (j in 1:S){
    d.mat[i, j] <- sqrt((ll[i, "Lat"] - ll[j, "Lat"])^2 + (ll[i, "Long"] - ll[j, "Long"])^2)*108
    # Q for self: where does the 108 come from there? Something about lat-long and distance in KM in ethiopia?
  }
}
rownames(d.mat) <- ll[, 1]
colnames(d.mat) <- ll[, 1]
d.mat <- d.mat/100





### 3. El Nino data ###

# Read in Nino 3.4 index:
# KV - this index may come from here: http://www.cgd.ucar.edu/cas/catalog/climind/TNI_N34/
# data from 1871 to present, V1 is the year
nino <- read.table("nino34.long.data", colClasses=rep("numeric", 13))

# vectorize data into one column
nino <- matrix(t(as.matrix(nino[, -1])), ncol=1)

# Make all dates, 12 months from 1871 to 2010
nino.dates <- paste(rep(1871:2010, each=12), unique(months(date.string)), sep=" ")

# cut dates and data off at April 2010
# nino only goes through April 2010 -- the last 8 values were "-99.99"
nino <- nino[1:which(nino.dates == "2010 April")]
nino.dates <- nino.dates[1:which(nino.dates == "2010 April")]

# subselect just the months that align with the rainfall data
# Start 3 months before first rainfall month
# We'll use nino with a 3-month lag
# Check in this, but in EDA I think the highest correlation with rainfall was
# to use nino from 3 months prior.
w <- which(nino.dates == "1991 October")
nino <- nino[w:length(nino.dates)]
nino.dates <- nino.dates[w:length(nino.dates)]

# center nino at its mean across these 223 months:
nino <- nino - mean(nino)
month.mat <- matrix(0, N, 12)
# repeat month data 18 times, add in 6 more months
# this: c(rep(month.days, 18), month.days[1:6], 28), constructs 18 repeats of the months, 
# plus the first 6 months plus one more element of 28
# so nino[1:223] are the values for each of 223 months, which is the latter expression, 
# so each value of nino will be repeated the number of days in the 223 months

# number of days per month for a non-leap year:
month.days <- c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)

# Align Nino start date with 3 month lag to data's start date:
month.mat[cbind(1:N, month.vec)] <- rep(nino, c(rep(month.days, 18), month.days[1:6], 28))
# month.mat contains the 3-month-lagged, normalized nino temperatures for each day in the data





### 4. ARC information ###

# Set up ARC indicator variable for each series:
arc <- c(1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0)

# ARC indicator (X^\text{ARC} in the paper)
# KV-ARC data is http://journals.ametsoc.org/doi/abs/10.1175/JCLI-D-12-00206.1
# KV-"(ARC) project aims to create an independent climate data record of sea surface 
# temperatures (SSTs) covering recent decades that can be used for climate change analysis"
X.arc <- as.list(rep(NA, S))
for (s in 1:S){
  X.arc[[s]] <- matrix(0, J[s], 1)
  X.arc[[s]][1, 1] <- 1
}






### 5. Design matrix for regression ###

# Create periodic predictor variables:
# KV - for each column (10 all together) he is creating a sine 
# and a cosine
# In one year, the wave goes up and down once, plot(sin.mat[,1])
M <- 10
m.seq <- c(1:10)
sin.mat <- matrix(NA, N, M)
cos.mat <- matrix(NA, N, M)
for (m in 1:M){
  sin.mat[, m] <- sin(2*pi*(1:N)*m.seq[m]/365)
  cos.mat[, m] <- cos(2*pi*(1:N)*m.seq[m]/365)
}

# Create X_t, which is the same across sites:
# KV - might be building a regression matrix
X <- cbind(rep(1, N), (1:N)/N - 0.5, ((1:N)/N - 0.5)^2, sin.mat[, 1:4], cos.mat[, 1:4], month.mat)

# KV-P is the number of trends in the time trend variable
P <- dim(X)[2]
X.names <- c("Intercept", "Time (linear)", "Time (quadratic)", "sin-1", "sin-2", "sin-3", "sin-4",
             "cos-1", "cos-2", "cos-3", "cos-4", paste(month.names, "Nino", sep="-"))






### 6. index the NA elements of the data set ###

# KV - he is creating a matrix of ones and zeros somehow related to 
# the NAs in the original data set, data
# KV - so na.mat is a list of 6 matrices ... the first 5 matrices will be 
# (6779 x 2) and the 6th will be (6779 x 5) based on groupings of the site data

# membership of each series within locations:
site.mat <- cbind(c(1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 6, 6, 6),
                  c(1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 3, 4, 5))

na.mat <- as.list(rep(NA, S))  # 1 = missing
for (s in 1:S){
  na.mat[[s]] <- matrix(0, J[s], N)
  na.mat[[s]][is.na(data[site.mat[, 1] == s, ])] <- 1
}


# KV: Some later parameters that are useful to save
#K <-3

#G<- 2000

#adapt <- 500

#burn <- 1001:2000

# Save the necessary data objects in a named list:
input.data <- list(Y = Y,  # observed rainfall
                   site.names = site.names,  # names of sites
                   N = N,  # number of time points
                   X.arc = X.arc,  # list of ARC indicators
                   S = S,  # number of locations
                   J = J,  # number of time series per location
                   date.string = date.string,  # vector of dates
                   month.vec = month.vec,  # vector of month indices
                   year.vec = year.vec,  # vector of year indices
                   d.mat = d.mat,  # distance matrix
                   X = X,  # design matrix for spatial regression
                   month.names = month.names, #month names for graphs
                   X.names = X.names,
                   P = P,  # number of predictors in X
                   na.mat = na.mat  # list of NA indicators for observed rainfall
                   )

save(input.data, file="input_data.RData")



