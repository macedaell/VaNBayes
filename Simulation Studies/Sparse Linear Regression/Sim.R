rm(list=ls())

setwd("S:\\Desktop\\SparseLinModel\\Sim\\")
source("adam.R")
source("hetnorm2.R")
source("logistic_regression2.R")
source("SSVS.R")
 
# These are the summary statistics, Z
summary_stats <- function(Y,X,H){
  b   <- as.vector(H%*%Y)
  sig <- sd(Y-X%*%b)
return(c(b/sig,log(sig),log(sd(b))))}

################      SIMULATION SETTINGS     ################ 
library(MASS)

nsims      <- 200      # number of simulated datasets
type       <- 1        # Type (1-6) of network architecture 
n          <- 50       # Sample size of each dataset
p          <- 20       # Number of covariates
beta0      <- rep(0,p) # True value of the parameters
beta0[1:2] <- 0.5
beta0[6]   <- 0.5
alpha0     <- 0
sigma0     <- 1
rho        <- 0.5      # Correlation between covariates
tau        <- 1        # Priors (see latex)
a          <- 0.5
b          <- 0.05
c          <- 2
d          <- 2
burn       <- 10000  # MCMC burn in
iters      <- 50000  # MCMC iterations 
m          <- 100000 # Number of synthetic datasets for training the DNN

if(type==1){nodes <- c( 50,10)}
if(type==2){nodes <- c( 50,25)}
if(type==3){nodes <- c(100,10)}
if(type==4){nodes <- c(100,25)}
if(type==5){nodes <- c(200,10)}
if(type==6){nodes <- c(200,25)}
print(c(p,type))

# Generate X (training) and Xp (testing) once

 # X for the real dataset
 set.seed(919)
 X          <- matrix(rnorm(n*p),n,p)
 for(j in 2:p){X[,j] <- rho*(X[,j-1]) + sqrt(1-rho^2)*X[,j]}
 X1         <- cbind(1,X)  
 H          <- ginv(t(X1)%*%X1)%*%t(X1)

 # X for the prediction set
 set.seed(0820)
 np         <- 10
 Xp         <- matrix(rnorm(np*p),np,p)
 for(j in 2:p){Xp[,j] <- rho*(Xp[,j-1]) + sqrt(1-rho^2)*Xp[,j]}

# Train our method

 set.seed(1)
 p_in  <- rbeta(m,c,d)
 sigma <- 1/sqrt(qgamma(runif(m,0.001,0.999),a,b)) # Prevent outliers
 alpha <- rnorm(m,0,tau)
 Z     <- matrix(0,m,p+3)
 beta  <- matrix(0,m,p)
 Yp    <- matrix(0,m,np)
 W     <- list()
 for(s in 1:m){
   set.seed(s)
   B        <- rbinom(p,1,p_in[s])*rnorm(p,0,tau)
   Ys       <- rnorm(n,alpha[s]+X%*%B,sigma[s])   
   Yp[s,]   <- rnorm(np,alpha[s]+Xp%*%B,sigma[s])   
   beta[s,] <- B
   Z[s,]    <- summary_stats(Ys,X1,H)
 }

 # Scale the inputs to [-1,1]
  qlev    <- seq(0,1,length=1001)
  Q       <- apply(Z,2,quantile,qlev)
  Z       <- quantile_scale(Z,Q)

 # Estimate the logistic regression weights 
  W_log <- list()
  for(j in 1:p){
    print(paste("Training covariate",j))
    gammaj     <- ifelse(beta[,j]==0,0,1)
    initj      <- init_logistic(p=ncol(Z),L1=nodes[1],L2=nodes[2],
                                init_mn=logit(mean(gammaj)))
    modelj     <- adam(w=initj, x=Z, y=gammaj, 
                       loss=loss_logistic, grad=grad_logistic)
    W_log[[j]] <- as.list(modelj$w)
  }

 # Estimate the heterogeneous normal weights 
  gamma  <- log(sigma)
  init   <- init_hetnorm(ncol(Z),nodes[1],nodes[2],
                         init_mn=mean(gamma),init_sd=sd(gamma))
  W_sig  <- adam(w=init, x=Z, y=gamma, 
                 loss=loss_hetnorm, grad=grad_hetnorm)$w 

 # Estimate the prediction weights 
  W_pred <- list()
  for(j in 1:np){
    print(paste("Training prediction",j))
    gammaj      <- Yp[,j]
    initj       <- init_hetnorm(p=ncol(Z),L1=nodes[1],L2=nodes[2],
                                init_mn=mean(gammaj),init_sd=sd(gammaj))
    modelj      <- adam(w=initj, x=Z, y=gammaj, 
                        loss=loss_hetnorm, grad=grad_hetnorm)
    W_pred[[j]] <- as.list(modelj$w)
  }


#Evaluate fits

 set.seed(2)
 p_in_test  <- rbeta(m,c,d)
 sigma_test <- 1/sqrt(qgamma(runif(m,0.001,0.999),a,b))
 alpha_test <- rnorm(m,0,tau)
 Z_test     <- matrix(0,m,ncol(Z))
 beta_test  <- matrix(0,m,p)
 Yp_test    <- matrix(0,m,np)
 for(s in 1:m){
   B             <- rbinom(p,1,p_in_test[s])*rnorm(p,0,tau)
   Ys            <- rnorm(n ,alpha_test[s]+X%*%B,sigma_test[s])   
   Yp_test[s,]   <- rnorm(np,alpha_test[s]+Xp%*%B,sigma_test[s])   
   beta_test[s,] <- B
   Z_test[s,]    <- summary_stats(Ys,X1,H)
 }
 Z_test <- quantile_scale(Z_test,Q)

 # PIPs
 CE <- CA <- BS <- rep(0,p)
 for(j in 1:p){
   gammaj <- ifelse(beta_test[,j]==0,0,1)
   pred   <- predict_logistic(W_log[[j]],Z_test)
   CE[j]  <- -mean(dbinom(gammaj,1,pred,log=TRUE))
   CA[j]  <- mean(gammaj == ifelse(pred>0.5,1,0))
   BS[j]  <- mean((gammaj-pred)^2)
 }
 
 # sigma
 pred   <- predict_hetnorm(W_sig,Z_test)
 LS_sig <- mean(dnorm(log(sigma_test),pred$mu,pred$sigma,log=TRUE))

 # Predictions
 LS_pred <- 0
 for(j in 1:np){
   pred    <- predict_hetnorm(W_pred[[j]],Z_test)
   LS_pred <- LS_pred + mean(dnorm(Yp_test[,j],pred$mu,pred$sigma,log=TRUE))/np
 }
 

# Start the simulation
 q_eval  <- c(0.05,0.50,0.95)
 pip0    <- pip1    <- matrix(0,nsims,p)
 q_sig0  <- q_sig1  <- matrix(0,nsims,length(q_eval))
 q_pred0 <- q_pred1 <- array(0,c(nsims,np,length(q_eval))) 
 Y_p     <- matrix(0,nsims,np)
 
for(sim in 1:nsims){

  set.seed(sim*9190)
  print(paste0("Starting dataset #",sim))

  # Generate data
   Y         <- rnorm(n,alpha0+X%*%beta0,sigma0)
   Y_p[sim,] <- rnorm(np,alpha0+Xp%*%beta0,sigma0)
   Z0        <- summary_stats(Y,X1,H)
   for(k in 1:ncol(Z)){Z0[k] <- stdq(Z0[k],Q[,k])}
 
  # Fit MCMC  
   fit0           <- SSVS(Y,X,xp=Xp,tau=tau,a=a,b=b,c=c,d=d,burn=burn,iters=iters,
                          init_beta=beta0,init_int=alpha0,init_sigma2=sigma0^2)
   pip0[sim,]     <- colMeans(fit0$beta[burn:iters,]!=0)
   q_sig0[sim,]   <- quantile(fit0$sigma[burn:iters],q_eval)
   q_pred0[sim,,] <- t(apply(fit0$pred[burn:iters,],2,quantile,q_eval))

  # Apply our method
   for(j in 1:p){pip1[sim,j] <- predict_logistic(W_log[[j]],Z0)}
   pred         <- predict_hetnorm(W_sig,Z0)
   q_sig1[sim,] <- exp(qnorm(q_eval,pred$mu,pred$sigma))   
   for(j in 1:np){
     pred            <- predict_hetnorm(W_pred[[j]],Z0)
     q_pred1[sim,j,] <- qnorm(q_eval,pred$mu,pred$sigma)
   }

  # Compare PIPs
   plot(pip0[sim,],pip1[sim,],xlim=0:1,ylim=0:1,main=sim)
   abline(0,1)
}

# Save output
 rm(Z,Z_test,fit0,Yp,Yp_test)
 save.image(paste0("Sim_p",p,"fit",type,".RData"))

# Summarize the results

p0 <- as.vector(pip0) 
p1 <- as.vector(pip1) 
plot(p0,p1,xlim=0:1,ylim=0:1,xlab="PIP MCMC",ylab="PIP NN")
abline(0,1)

print("p, nodes")
print(p)
print(nodes)
print("CE, CA, BS for PIP")
print(mean(CE))
print(mean(CA))
print(mean(BS))
print("LS for sigma")
print(LS_sig)
print("LS for prediction")
print(LS_pred)
print("Cor PIP, MCMC/DNN")
print(cor(p0,p1))
print("Cor median sigma, MCMC/DNN")
print(cor(q_sig1[,2],q_sig0[,2]))
print("Cor median yp, MCMC/DNN")
print(cor(as.vector(q_pred1[,,2]),as.vector(q_pred0[,,2])))
