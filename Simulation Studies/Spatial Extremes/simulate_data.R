rm(list=ls())
library(MASS)
library(spam)
library(SpatialExtremes)
library(evd)

#setwd("S:\\Desktop\\MSP\\")
source("adam.R")
source("hetnorm2.R")
source("simu_Dombry_et_al.R")

orig2log = function(data, param){
  if(param =="nu"){
    out = log(data/(2 - data))
  }
  if(param =="range"){
    out = log(data)
  }
  return(out)
}
log2orig = function(data, param){
  if(param =="nu"){
    out = 2*exp(data)/(1+ exp(data))
  }
  if(param =="range"){
    out = exp(data)
  }
  return(out)
}

################      SIMULATION SETTINGS     ################
n1         <- 100   # Sample size of each dataset is n1^2
n.size <-  10       # test data grid length
m          <- 10000 # Number of synthetic datasets for training the DNN
sd.prior   <- 1

n.rep  = 50
n.bins = 10
np     <-  5

x.test <- runif(n1, 0, n.size)
y.test <- runif(n1, 0, n.size)
coord  <- cbind(x.test, y.test)
save(coord, file="data/coords_100loc_size10.Rdata")

set.seed(1)

theta                 <- matrix(NA,m,np)
theta[,c(1,2,4)]      <- matrix(rnorm(m*3, 0, sd.prior),m,3)
theta[,3]      <- rnorm(m, 1, 0.5)
theta[,5]      <- rnorm(m, 0, 0.1)

theta.orig <- cbind(log2orig(param = "range", data=theta[,1]),
                    log2orig(param =    "nu", data=theta[,2]),
                    theta[,3],
                    log2orig(param = "range", data=theta[,4]),
                    theta[,5])

Y = array(NA, c(n.rep, nrow(coord), m))
Z = matrix(NA, m, n.bins+6+3)

for(s in 1:m){

  vario <- function(x) ((1/theta.orig[s,1])*(sqrt(sum(x^2))))^theta.orig[s,2]

  if(theta.orig[s,1]>0.1){  # Can we remove this?

    Y[,,s]     <- simu_specfcts(model="brownresnick",
                                loc=theta.orig[s,3], scale=theta.orig[s,4], shape=theta.orig[s,5],
                                no.simu=n.rep,
                                coord=coord,
                                vario=vario)$res

    Z[s,1:n.bins] <- fmadogram(Y[,,s], coord, n.bins=n.bins)[1:n.bins,3]
    Z[s,(n.bins+1):(n.bins+6)] <-  quantile(c(Y[,,s]), c(0.5, 0.7, 0.9, 0.95, 0.99, 1))
    Z[s,(n.bins+7):(n.bins+9)] <-  fgev(c(Y[,,s]))$estimate[1:3] # xi, mu, beta

  }
  if(s%%100==0){print(s)}
}

junk       <- is.na(rowSums(Z))
theta      <- theta[!junk,]
theta.orig <- theta.orig[!junk,]
Z          <- Z[!junk,]
Y          <- Y[,,!junk]

save(Y, file="data/Y_gev_thetanorm01all5_sumall.Rdata")
save(Z, file="data/Z_gev_thetanorm01all5_sumall.Rdata")
save(theta, file="data/theta_gev_thetanorm01all5_sumall.Rdata")
