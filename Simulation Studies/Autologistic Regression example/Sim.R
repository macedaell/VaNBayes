rm(list=ls())
library(MASS)
library(spam)

# This is the main function that runs the simulation study
setwd("S:\\Desktop\\Autologistic\\Sim\\")
source("adam.R")
source("hetnorm2.R")
source("autologistic.R")
 
# These are the summary statistics, Z
summary_stats <- function(Y,X,N1,N2,N3){
    bhat  <- glm(Y~X-1,family="binomial")$coef
    res  <- Y-expit(X%*%bhat)
    Z    <- c(bhat,
              log(GearyC(res,N1)),
              log(GearyC(res,N2)),
              log(GearyC(res,N3)))
return(Z)}

################      SIMULATION SETTINGS     ################ 

nsims      <- 1000     # number of simulated datasets
type       <- 6        # Type (1-4) of network architecture 
g_type     <- 2        # Level (1-2) of gamma
n1         <- 20       # Sample size of each dataset is n1^2
p          <- 5        # Number of covariates
np         <- p+1
nz         <- p+3
beta0      <- rep(0,p) # True value of the parameters
beta0[2]   <- 0.5
gamma0     <- ifelse(g_type==1,0.25,0.50)
theta0     <- c(beta0,log(gamma0))
m          <- 100000 # Number of synthetic datasets for training the DNN

if(type==1){nodes <- c( 50,10)}
if(type==2){nodes <- c( 50,25)}
if(type==3){nodes <- c(100,10)}
if(type==4){nodes <- c(100,25)}
if(type==5){nodes <- c(200,10)}
if(type==6){nodes <- c(200,25)}
print(c(p,type))

# Set up the spatial network
n <- n1^2
s <- expand.grid(1:n1,1:n1)
d <- as.matrix(dist(s))
N1 <- N2 <- N3 <- list()
for(i in 1:n){
  N1[[i]] <- which(d[i,]==1)
  N2[[i]] <- which(d[i,]==sqrt(2))
  N3[[i]] <- which(d[i,]==2)
}
N     <- as.spam(ifelse(d==1,1,0))
block <- ((s[,1]-s[,2])%%2==0) + 1

# Generate X

 # X for the real dataset
 set.seed(919)
 X      <- matrix(rnorm(n*p),n,p)
 X[,1]  <- 1

# Train our method
  print("Generating training data")
  set.seed(1)
  theta <- matrix(rnorm(m*np),m,np)
  Z     <- matrix(NA,m,nz)
  W     <- list()
  for(s in 1:m){
   set.seed(s)
   kappa    <- expit(X%*%theta[s,1:p])
   gamma    <- exp(theta[s,p+1])
   Ys       <- rautologistic(kappa,gamma,N1)
   if(sum(Ys) > 10 & sum(1-Ys)>10){
      Z[s,]    <- summary_stats(Ys,X,N1,N2,N3)
   }
  }
  junk  <- is.na(Z[,1])
  theta <- theta[!junk,]
  Z     <- Z[!junk,]

 # Scale the inputs to [-1,1]
  qlev    <- seq(0,1,length=1001)
  Q       <- apply(Z,2,quantile,qlev)
  Z       <- quantile_scale(Z,Q)

 # Estimate the weights 
  W  <- list()
  for(j in 1:np){
    print(paste("Training parameter",j))
    gammaj     <- theta[,j]
    initj      <- init_hetnorm(p=ncol(Z),L1=nodes[1],L2=nodes[2],
                               init_mn=mean(gammaj),init_sd=sd(gammaj))
    modelj     <- adam(w=initj, x=Z, y=gammaj, 
                       loss=loss_hetnorm, grad=grad_hetnorm)
    W[[j]] <- as.list(modelj$w)
  }

#Evaluate fits

 print("Generating testing data")
 set.seed(2)
 theta_test  <- matrix(rnorm(m*np),m,np)
 Z_test      <- matrix(NA,m,nz)
 for(s in 1:m){
   kappa      <- expit(X%*%theta_test[s,1:p])
   gamma      <- exp(theta_test[s,p+1])
   Ys         <- rautologistic(kappa,gamma,N1)
   if(sum(Ys) > 10 & sum(1-Ys)>10){
     Z_test[s,] <- summary_stats(Ys,X,N1,N2,N3)
   }
 }
 junk       <- is.na(Z_test[,1])
 theta_test <- theta_test[!junk,]
 Z_test     <- Z_test[!junk,]
 Z_test     <- quantile_scale(Z_test,Q)

 # Fits
 LS  <- rep(0,np)
 PIT <- matrix(0,nrow(Z_test),np)
 for(j in 1:np){
   pred    <- predict_hetnorm(W[[j]],Z_test)
   LS[j]   <- mean(dnorm(theta_test[,j],pred$mu,pred$sigma,log=TRUE))
   PIT[,j] <- pnorm(theta_test[,j],pred$mu,pred$sigma)
 }
 plot(NA,xlim=0:1,ylim=0:1,xlab="Theoretical quantile",ylab="Observed quantile")
 for(j in 1:np){
    qqq   <- seq(0.01,0.99,0.01)
    lines(qqq,quantile(PIT[,j],qqq),lty=j)
 }

# Start the simulation
 q_eval  <- c(0.05,0.50,0.95)
 q       <- q_pred1 <- array(0,c(nsims,np,length(q_eval))) 
 
for(sim in 1:nsims){

  set.seed(sim*9190)
  print(paste0("Starting dataset #",sim))

  # Generate data
   kappa  <- expit(X%*%beta0)
   Y      <- rautologistic(kappa,gamma0,N1)
   Z0     <- summary_stats(Y,X,N1,N2,N3)
   for(k in 1:nz){Z0[k] <- stdq(Z0[k],Q[,k])}

  # Apply our method
   for(j in 1:np){
     pred      <- predict_hetnorm(W[[j]],Z0)
     q[sim,j,] <- qnorm(q_eval,pred$mu,pred$sigma)
   }
   matplot(q[sim,,],col=1,lty=c(2,1,2),type="l",xlab="Parameter",ylab="Posterior",main=sim)
   lines(theta0,col=2,lwd=2)
}

# Save output
 rm(Z,Z_test)
 save.image(paste0("Sim_g",g_type,"fit",type,".RData"))

# Summarize the results
cov <- rep(0,np)
for(j in 1:np){
  cov[j] <- mean((q[,j,1]<theta0[j]) & (theta0[j]<q[,j,3]))
} 
save.image(paste0("Sim_g",g_type,"fit",type,".RData"))

print("p, nodes")
print(p)
print(nodes)
print("LS")
print(LS)
print("Coverage")
print(round(cov,2))

boxplot(q[,,2],xlab="Parameter",ylab="Post median")
   lines(theta0,col=2,lwd=2)
