
# Elliot Maceda
# NCSU, May 2024


###################### PACKAGES AND FUNCTIONS REQUIRED #########################
library(ggplot2)
library(MASS)
library(spam)

# variational posterior
source("hetnorm2.R")

# neural network optimizer
source("adam.R")


########## IMPORTING THE TRAINING DATA FROM THE TRAINING DISTRIBUTION ##########

imported_data = as.matrix(read.csv("training_data.csv"))[,-1]
Y = imported_data[,-c(1,2,3,4)]
theta = imported_data[,c(1,2,3,4)]


##################### IMPORTING THE OBSERVED DATA ##############################

# importing this for the log-density of each population
load("zika Brazil data.rdata")

true_data = as.matrix(read.csv("boserved_data.csv"))[,-1]

true_data_mat = matrix(true_data, ncol = 40)

########################### SUMMARY STATISTICS #################################

# Neural Network Architecture -- completely connected with two layers
neural_network_size = "small"
summary_vars = 4 # 4 (83%), 6 (90%), 9 (95%), 20 (99%), 30 (99.5%)

if (neural_network_size == "small"){
  layer1_nodes = 50
  layer2_nodes = 10
} else if (neural_network_size == "medium"){
  layer1_nodes = 100
  layer2_nodes = 25
} else if (neural_network_size == "large"){
  layer1_nodes = 150
  layer2_nodes = 50
} else {
  layer1_nodes = 200
  layer2_nodes = 100
}

# # Use PCA to summarize the Y dataset
number_of_pcs = summary_vars
eigen_decomp = eigen(cov(Y))
Z = Y%*%eigen_decomp$vectors[,1:number_of_pcs]
true_Z = true_data%*%eigen_decomp$vectors[,1:number_of_pcs]

reconstructed_true_data = true_Z%*%t(eigen_decomp$vectors[,1:number_of_pcs])
reconstructed_true_data_mat = matrix(reconstructed_true_data, ncol = 40)

# standardize the PC scores for the neural network
Z_combined = rbind(true_Z, Z)
std_Z_combined = scale(Z_combined)
std_true_Z = std_Z_combined[1,]
std_Z = std_Z_combined[-1,]

######################### FITTING THE NEURAL NETWORK ###########################

# endpoints of the uniform distributions
beta0_lower = -3
beta0_upper = 1
beta1_lower = -1
beta1_upper = 1
nu_lower = 1.01
nu_upper = 10


# transform the parameters so the heterogeneous normal model makes sense for them
normal_theta = theta
normal_theta[,1] = qnorm((theta[,1] - (beta0_lower))/(beta0_upper-(beta0_lower))) # beta0
normal_theta[,2] = qnorm((theta[,2] - (beta1_lower))/(beta1_upper-(beta1_lower))) # beta1
normal_theta[,3] = log(theta[,3])                                                 # phi
normal_theta[,4] = qnorm((theta[,4] - (nu_lower))/(nu_upper-(nu_lower)))          # nu


# fit the neural networks (one for each parameter)
W  <- list()
for(j in 1:ncol(theta)){
  print(paste("Training parameter",j))
  
  gammaj     <- normal_theta[,j]# theta[,j]# normal_theta[,j] # response for the jth neural network
  

  initj      <- init_hetnorm(p=ncol(std_Z),L1=layer1_nodes,L2=layer2_nodes,
                             init_mn=mean(gammaj),init_sd=sd(gammaj)) # initialize weights with normal distribution
  
  # adaptive SGD using our standardized PC Scores and the log-normal as the loss function
  modelj     <- adam(w=initj, x=std_Z, y=gammaj, priorweights = 1, lr = 1e-2, early_stop = 20,
                     loss=loss_hetnorm, grad=grad_hetnorm, verbose = 1, batchsize = 100) 
  
  W[[j]] <- as.list(modelj$w) # save the weights of each neural network
}


########################### Estimating the Posterior ###########################


# Fits
# NOTE: Here, we use std_Z and normal_theta
# however, in reality these need to be generated independently of the training data

# log-scores
LS  <- rep(0,4)
PIT <- matrix(0,nrow(std_Z),4)
for(j in 1:4){
  pred    <- predict_hetnorm(W[[j]],std_Z)
  LS[j]   <- mean(dnorm(normal_theta[,j],pred$mu,pred$sigma,log=TRUE))
  PIT[,j] <- pnorm(normal_theta[,j],pred$mu,pred$sigma)
}

# PIT plots
plot(NA,xlim=0:1,ylim=0:1,xlab="Theoretical quantile",ylab="Observed quantile")
mycolors = c("red", "blue", "green")
for(j in 1:3){
  qqq   <- seq(0.01,0.99,0.01)
  lines(qqq,quantile(PIT[,j],qqq),lty=j, col = mycolors[j])
}
legend(x = c(0.8,1.00), y= c(0,0.2), legend = c(expression(beta[0]), expression(beta[1]), expression(phi)), col = mycolors, lty = c(1,2,3))


# estimating the variational posteriors corresponding to our data
pred      <- predict_hetnorm(W[[1]],std_true_Z)
normal_posterior = rnorm(10000,mean = pred$mu, sd = pred$sigma)
beta0_posterior = pnorm(normal_posterior)*(beta0_upper-(beta0_lower)) + (beta0_lower)

pred      <- predict_hetnorm(W[[2]],std_true_Z)
pred$mu
normal_posterior = rnorm(10000,mean = pred$mu, sd = pred$sigma)
beta1_posterior = pnorm(normal_posterior)*(beta1_upper-(beta1_lower)) + (beta1_lower)

pred      <- predict_hetnorm(W[[3]],std_true_Z)
pred$mu
normal_posterior = rnorm(10000,mean = pred$mu, sd = pred$sigma)
phi_posterior = exp(normal_posterior)

pred      <- predict_hetnorm(W[[4]],std_true_Z)
pred$mu
normal_posterior = rnorm(10000,mean = pred$mu, sd = pred$sigma)
nu_posterior = pnorm(normal_posterior)*(nu_upper-(nu_lower)) + (nu_lower)

# variational posterior means
print(mean(beta0_posterior))
print(mean(beta1_posterior))
print(mean(phi_posterior))

# 95% credible intervals
print(c(quantile(beta0_posterior, 0.025), quantile(beta0_posterior, 0.975)))
print(c(quantile(beta1_posterior, 0.025), quantile(beta1_posterior, 0.975)))
print(c(quantile(phi_posterior, 0.025), quantile(phi_posterior, 0.975)))

############################ PCA IMAGE PICS ####################################

true_data_mat_ordered = true_data_mat[order(log_density, decreasing = T), ]
reconstructed_true_data_mat_ordered = reconstructed_true_data_mat[order(log_density, decreasing = T), ]

true_data_mat_rotated <- apply(true_data_mat_ordered, 2, rev)
image(1:ncol(true_data_mat_ordered), 1:nrow(true_data_mat_ordered), t(true_data_mat_rotated), ylab = "State (from most to least dense)", xlab = "Week", yaxt = "n", col = hcl.colors(12, "viridis", rev = FALSE), zlim = c(min(reconstructed_true_data_mat_ordered, max(true_data_mat_ordered))))

reconstructed_true_data_mat_rotated <- apply(reconstructed_true_data_mat_ordered, 2, rev)
image(1:ncol(reconstructed_true_data_mat_ordered), 1:nrow(reconstructed_true_data_mat_ordered), t(reconstructed_true_data_mat_rotated), ylab = "State (from most to least dense)", xlab = "Week", yaxt = "n", col = hcl.colors(12, "viridis", rev = FALSE), zlim = c(min(reconstructed_true_data_mat_ordered, max(true_data_mat_ordered))))

loading_matrix = eigen_decomp$vectors[,1:number_of_pcs]

for (k in 1:number_of_pcs){
  loading_weights = loading_matrix[,k]
  loading_weights_mat = matrix(loading_weights, ncol = 40, byrow = F)
  loading_weights_mat_ordered = loading_weights_mat[order(log_density, decreasing = T), ]
  
  loading_weights_mat_ordered
  true_data_mat_rotated <- apply(loading_weights_mat_ordered, 2, rev)
  image(1:ncol(loading_weights_mat_ordered), 1:nrow(loading_weights_mat_ordered), t(true_data_mat_rotated), ylab = "State", xlab = "Week", yaxt = "n", col = hcl.colors(12, "viridis", rev = FALSE))
}


