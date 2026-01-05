rm(list=ls())
library(MASS)
library(spam)
library(SpatialExtremes)
library(evd)

# This is the main function that runs the simulation study

source("adam.R")
source("hetnorm2.R")
source("simu_Dombry_et_al.R")

load("data/Y_gev_thetanorm01all5_sumall.Rdata")
load("data/Z_gev_thetanorm01all5_sumall.Rdata")
load("data/theta_gev_thetanorm01all5_sumall.Rdata")
load("data/coords_100loc_size10.Rdata")

m          <- 10000 # Number of synthetic datasets for training the DNN
sd.prior   <- 1
np     <-  5

mnZ = apply(Z, 2, mean)
sdZ = apply(Z, 2, sd)

for(j in 1:ncol(Z)){Z[,j] <- pnorm(Z[,j],mnZ[j],sdZ[j])}

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

type=4
if(type==1){nodes <- c( 50,10)}
if(type==2){nodes <- c( 50,25)}
if(type==3){nodes <- c(100,10)}
if(type==4){nodes <- c(100,25)}
if(type==5){nodes <- c(200,10)}
if(type==6){nodes <- c(200,25)}

# Estimate the weights
W  <- list()
for(j in 1:np){
  print(paste("Training parameter",j))
  gammaj <- theta[,j]
  initj  <- init_hetnorm(p=ncol(Z[,1:16]),L1=nodes[1],L2=nodes[1],
                         init_mn=mean(gammaj),init_sd=sd(gammaj))
  modelj <- adam(w=initj, x=Z[,1:16], y=gammaj,
                 loss=loss_hetnorm, grad=grad_hetnorm)
  W[[j]] <- as.list(modelj$w)
}

# set variables
np     <-  5

theta01.orig  <- c(1.0, 2)
theta02.orig  <- c(1.0, 1.5)
theta03.orig  <- c(1.0, 1.5)
theta04.orig  <- c(1.0, 2)
theta05.orig  <- c(-0.1, 0.1)

theta0.orig = as.matrix(expand.grid(theta01.orig, theta02.orig,
                                    theta03.orig, theta04.orig, theta05.orig))

ntheta      <- nrow(theta0.orig)

theta0 <- theta0.orig
for(j in 1:ntheta){
  theta0[j,] = c(orig2log(param = "range", data=theta0.orig[j,1]),
                 orig2log(param =    "nu", data=theta0.orig[j,2]),
                          theta0.orig[j,3],
                 orig2log(param = "range", data=theta0.orig[j,4]),
                          theta0.orig[j,5])
}

# Start the simulation
nsims      <- 50   # number of simulated datasets
n.rep  = 50
n.bins = 10

q_eval  <- c(0.05,0.50,0.95)
q       <- array(0,c(nsims,ntheta,np,length(q_eval)))

for(k in 1:ntheta){for(sim in 1:nsims){
  set.seed(sim*9190)
  print(paste0("Starting dataset #",sim))

  # Generate data
  vario <- function(x) ((1/theta0.orig[k,1])*(sqrt(sum(x^2))))^theta0.orig[k,2]
  Y0     <- simu_specfcts(model="brownresnick",
                          loc=theta0.orig[k,3], scale=theta0.orig[k,4], shape=theta0.orig[k,5],
                          no.simu=n.rep,
                          coord=coord,
                          vario=vario)$res
  Z0 = rep(NA, n.bins+6+2)
  Z0[1:n.bins] <- fmadogram(Y0, coord, n.bins=n.bins)[1:n.bins,3]
  Z0[(n.bins+1):(n.bins+6)] <-  quantile(c(Y0), c(0.5, 0.7, 0.9, 0.95, 0.99, 1))
  Z0[(n.bins+7):(n.bins+9)] <-  fgev(c(Y0))$estimate[1:3] # xi, mu, beta

  Z0    <- pnorm(Z0, mnZ,sdZ)

  # Apply our method
  for(j in 1:np){
    pred      <- predict_hetnorm(W[[j]],Z0[1:16])
    q[sim,k,j,] <- qnorm(q_eval,pred$mu,pred$sigma)
  }
  matplot(q[sim,k,,],col=1,lty=c(2,1,2),type="l",xlab="Parameter",ylab="Posterior",main=sim)
  lines(theta0[[k]],col=2,lwd=2)
}}

# Summarize the results
cov <- matrix(0,ntheta,np)
for(k in 1:ntheta){for(j in 1:np){
  cov[k,j] <- mean((q[,k,j,1]<theta0[k,j]) & (theta0[k,j]<q[,k,j,3]))
}}


png(file=paste0(save_dir, "/plots/boxplots_postmedian_thetanorm01all5_sumall_nomle2.png"),
    width=500, height=500)

par(mar=c(1, 4, 1, 1.2), mfrow=c(2,2), oma = c(1, 1, 0.5, 0.2))
ipc = c(25,5,9,29)
for(ip in ipc){
  labs <- c(expression(theta[1]), expression(theta[2]), expression(theta[3]),
            expression(theta[4]), expression(theta[5]))

  if(ip %in% c(25,9)){
    boxplot(q[,ip,,2],axes=F,outline=F, ylab="Posterior median",cex.lab=1.25, #main=ip,
            ylim = c(-0.6,1.9))
  } else{
    boxplot(q[,ip,,2],axes=F,outline=F, cex.lab=1.25, ylim = c(-0.6,1.9))
  }

  axis(2,cex.axis=1.25,cex.lab=1.25)
  axis(1,at=1:5,lab=labs,cex.axis=1.25,cex.lab=1.25)
  points(theta0[ip,],col=2,pch=19)
  text(1:5,-0.72,format(cov[ip,],nsmall = 2),pos=3,cex=1.25)
}

dev.off()


##############
## PIT histograms
################

# Start the simulation
z.len = 10000
m          <- 10000 # Number of synthetic datasets for training the DNN
p       <-matrix(NA, z.len, np)
sd.prior=1

# Generate data
for(k in 1:z.len){

  theta.pit = rep(NA, 5)
  theta.pit[c(1,2,4)] <- rnorm(3, 0, sd.prior)
  theta.pit[3]       <- rnorm(1, 1, 0.5)
  theta.pit[5]        <- rnorm(1, 0, 0.1)

  theta.pit.orig =  c(log2orig(param = "range", data=theta.pit[1]),
                 log2orig(param = "nu", data=theta.pit[2]),
                 theta.pit[3],
                 log2orig(param = "range", data=theta.pit[4]),
                 theta.pit[5])

  vario <- function(x) ((1/theta.pit.orig[1])*(sqrt(sum(x^2))))^theta.pit.orig[2]
  Y0     <- simu_specfcts(model="brownresnick",
                          loc=theta.pit.orig[3], scale=theta.pit.orig[4], shape=theta.pit.orig[5],
                          no.simu=n.rep,
                          coord=coord,
                          vario=vario)$res

  Z0 = rep(NA, n.bins+6+2)
  Z0[1:n.bins] <- fmadogram(Y0, coord, n.bins=n.bins)[1:n.bins,3]
  Z0[(n.bins+1):(n.bins+6)] <-  quantile(c(Y0), c(0.5, 0.7, 0.9, 0.95, 0.99, 1))
  Z0[(n.bins+7):(n.bins+9)] <-  fgev(c(Y0))$estimate[1:3] # xi, mu, beta

  Z0    <- pnorm(Z0, mnZ,sdZ)

  # Apply our method
  for(j in 1:np){
    pred      <- predict_hetnorm(W[[j]],Z0[,1:16])
    p[k,j] <- pnorm(theta.pit[j],pred$mu,pred$sigma)
  }
  if(k%%100==0){print(k)}
}

hist(p[,2])

qunif.eval = seq(0, 1, by=0.05)

q.pred.range = quantile(p[,1], probs = qunif.eval)
q.pred.smooth = quantile(p[,2], probs = qunif.eval)
q.pred.location = quantile(p[,3], probs = qunif.eval)
q.pred.sigma = quantile(p[,4], probs = qunif.eval)
q.pred.scale = quantile(p[,4], probs = qunif.eval)
q.unif = qunif(qunif.eval)

png(file=paste0(save_dir, "/plots/qqplot_gev_thetanorm01all5_sumall.png"), width=500, height=500)

plot(q.unif, q.pred.range,  type="l", xlab="Theoretical quantile", ylab="Observed quantile",
     cex.lab=1.4, cex.axis=1.4, lwd=2, lty=2, col=2)
points(q.unif,q.pred.smooth, type="l", col=3, lwd=2, lty=3)
points(q.unif,q.pred.location, type="l", col=4, lwd=2, lty=4)
points(q.unif,q.pred.sigma, type="l", col=5, lwd=2, lty=5)
points(q.unif,q.pred.scale, type="l", col=6, lwd=2, lty=6)
abline(0,1, lty=1, lwd=2)
legend("bottomright", c(expression(theta[1]),  expression(theta[2]),
       expression(theta[3]), expression(theta[4]), expression(theta[4])),
       col = 2:6, lwd=2, lty=2:6)

dev.off()

save.image(paste0("Sim_GEV",type,".RData"))
