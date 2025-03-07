
# Zika model fit via ABC
# use "zika Brazil.rdata" from main data analysis

rm(list = ls())
source("adam.R")
source("hetnorm2.R")
options(warn = 2)

#########
# Setup #
#########

# These are the summary statistics, Z
summary_stats <- function(Y,X,H){
  b   <- as.vector(H%*%Y)
  sig <- sd(Y-X%*%b)
  # return(c(log(sig), b, b/sig,log(sd(b))))}
  return(c(log(sig), b))}

################      SIMULATION SETTINGS     ################ 
library(MASS)

nsims      <- 100      # number of simulated datasets
n          <- 100       # Sample size of each dataset
p          <- 45       # would like this to be harder as we go up

# true beta values
beta0 = 2*((1:p %% 5) - 2)


rho        <- 0.0      #  high is bad for us
SNR = 0.8              # high is good for us
# Generate X (training) and Xp (testing) once
# X for the real dataset
set.seed(2025)
X          <- matrix(rnorm(n*p),n,p)
for(j in 2:p){X[,j] <- rho*(X[,j-1]) + sqrt(1-rho^2)*X[,j]}
X1         <- cbind(1,X)  
H          <- ginv(t(X1)%*%X1)%*%t(X1)
cat("X generated as", head(c(X)), "\n")
write.csv(X, file = paste0("./shared_data/p_",p,"/X_p",p,".csv"))

sigma0 <- as.numeric(sqrt(var(X%*%beta0)/SNR - var(X%*%beta0)))



alpha0     <- 1
a          <- 0.5 # used to be 0.5
b          <- 0.05
tau        <- .0005
SD         <- 10
N          <- 10000 # Number of synthetic datasets for training the DNN
type       <- 3        # Type (1-6) of network architecture 
my_lr      <- 0.0001
if(type==1){nodes <- c( 50,10)}
if(type==2){nodes <- c( 50,25)}
if(type==3){nodes <- c(100,10)}
if(type==4){nodes <- c(100,25)}
if(type==5){nodes <- c(200,10)}
if(type==6){nodes <- c(200,25)}
print(c(p,type))



# Draw from our Prior
sigma <- 1/sqrt(qgamma(runif(N,0.001,0.999),a,b)) # Prevent outliers
lambda <- sqrt(1/rbeta(N, 0.5, 0.5) - 1)
Z     <- matrix(0,N,p+2)
beta  <- matrix(0,N,p)
alpha  <- numeric(N)
Ys_kept    <- matrix(0,N,n)
for(s in 1:N){
  # set.seed(s)
  
  # generating_sd <- lambda[s]^2*tau^2 # - sparsity prior
  generating_sd <- SD # - wide normal
  
  B        <- rnorm(p,0,generating_sd)
  alpha[s]    <- rnorm(1,0,generating_sd)
  Ys       <- rnorm(n,alpha[s]+X%*%B,sigma[s])
  Ys_kept[s,] <- Ys
  beta[s,] <- B
  Z[s,]    <- summary_stats(Ys,X1,H)
}


# Scale the inputs to [-1,1]
qlev    <- seq(0,1,length=1001)
Q       <- apply(Z,2,quantile,qlev)
Z       <- quantile_scale(Z,Q)


# Estimate the regression weights 
W <- list()
for(j in 1:p){
  print(paste("Training covariate",j))
  gammaj     <- beta[,j]
  initj      <- init_hetnorm(p=ncol(Z),L1=nodes[1],L2=nodes[2],
                              init_mn=mean(gammaj), init_sd=sd(gammaj))
  modelj     <- adam(w=initj, x=Z, y=gammaj, lr = my_lr,
                     loss=loss_hetnorm, grad=grad_hetnorm)
  W[[j]] <- as.list(modelj$w)
}

# Estimate the intercept weights 
gamma  <- alpha
init   <- init_hetnorm(ncol(Z),nodes[1],nodes[2],
                       init_mn=mean(gamma),init_sd=sd(gamma))
W_alpha  <- adam(w=init, x=Z, y=gamma, lr = my_lr, 
               loss=loss_hetnorm, grad=grad_hetnorm)$w 


# Estimate the heterogeneous normal weights 
gamma  <- log(sigma)
init   <- init_hetnorm(ncol(Z),nodes[1],nodes[2],
                       init_mn=mean(gamma),init_sd=sd(gamma))
W_sig  <- adam(w=init, x=Z, y=gamma, lr = my_lr, 
               loss=loss_hetnorm, grad=grad_hetnorm)$w 

# save(W, W_alpha, W_sig, file = "successful_VaNBayes_fits.RData")

########################### SIMULATION STUDY ###################################

# Start the simulation
q_eval  <- c(0.05,0.50,0.95)
q_beta <- array(0,c(nsims,p,length(q_eval)))
q_alpha  <- matrix(0,nsims,length(q_eval))
q_sigma  <- matrix(0,nsims,length(q_eval))
sim_study_datasets     <- matrix(0,nsims,p+2)

for(sim in 1:nsims){
  
  # set.seed(sim*2025)
  print(paste0("Starting dataset #",sim))
  
  # Generate data
  Y   <- rnorm(n,alpha0+X%*%beta0,sigma0)
  # sim_study_datasets[sim,] <- Y
  Z0        <- summary_stats(Y,X1,H)
  
  sim_study_datasets[sim,] <- Z0
  for(k in 1:ncol(Z)){Z0[k] <- stdq(Z0[k],Q[,k])}
  
  # Apply our method
  for(j in 1:p){
    pred         <- predict_hetnorm(W[[j]],Z0)
    q_beta[sim,j,] <- qnorm(q_eval,pred$mu,pred$sigma)
  }
  pred         <- predict_hetnorm(W_alpha,Z0)
  q_alpha[sim,] <- qnorm(q_eval,pred$mu,pred$sigma) 
  pred         <- predict_hetnorm(W_sig,Z0)
  q_sigma[sim,] <- exp(qnorm(q_eval,pred$mu,pred$sigma))   
  
}
cat("Y generated as", head(c(sim_study_datasets)), "\n")
write.csv(sim_study_datasets, file = paste0("./shared_data/p_",p,"/sim_study_datasets_p",p,".csv"))


# coverage
beta_coverage = c()
for (j in 1:p){
  beta_coverage = c(beta_coverage, sum((q_beta[,j,1] < beta0[j]) & (beta0[j] < q_beta[,j,3]))/nsims)
}
alpha_coverage = sum((q_alpha[,1] < alpha0) & (alpha0 < q_alpha[,3]))/nsims
sigma_coverage = sum((q_sigma[,1] < sigma0) & (sigma0 < q_sigma[,3]))/nsims

# MAD
beta_MAD = c()
for (j in 1:p){
  beta_MAD = c(beta_MAD, mean(abs(q_beta[,j,2] - beta0[j])))
}
alpha_MAD = mean(abs(q_alpha[,2] - alpha0))
sigma_MAD = mean(abs(q_sigma[,2] - sigma0))


# boxplot of posterior medians vs. true values
labs <- 1:p
c1 <- beta_coverage
boxplot(q_beta[,,2],axes=F,outline=F, ylab="Posterior median",cex.lab=1.25)
axis(2,cex.axis=1.25,cex.lab=1.25)
axis(1,at=1:p,lab=labs,cex.axis=1.25,cex.lab=1.25)
points(beta0,col=2,pch=19)
text(1:p,-0.8,format(c1,nsmall=2),pos=3,cex=1.25)


beta_medians = q_beta[,,2]
alpha_medians = q_alpha[,2]
sigma_medians = q_sigma[,2]


save(beta_coverage, 
     alpha_coverage, 
     sigma_coverage, 
     beta_MAD, 
     alpha_MAD, 
     sigma_MAD,
     beta_medians,
     alpha_medians,
     sigma_medians,
     file = paste0("./results/VaNBayes_results_p",p,".RData"))


