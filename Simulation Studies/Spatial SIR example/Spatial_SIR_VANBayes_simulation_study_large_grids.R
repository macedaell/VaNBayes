
# Elliot Maceda
# NCSU, May 2024
# Spatial SIR Simulation Study

###################### PACKAGES AND FUNCTIONS REQUIRED #########################
library(ggplot2)
library(MASS)
library(spam)


setwd("C:/")

# Spatial SIR Simulator
source("Spatial_SIR_function_large_grids.R")

# variational posterior
source("autologistic example/hetnorm2.R")

# neural network optimizer
source("autologistic example/adam.R")

# creating the grid and variables that help run the Spatial SIR simulations efficiently
source("SIR_Grid_Setup.R")



########## Generating NN Training Data using Training Distribution #############

set.seed(NULL)
m <- 10000 # recommended 100000

# generate the parameters
betas = runif(m, 0.1, .9)
phis = runif(m, 0.1, .9)
etas = runif(m, 0.1, .9)
theta = cbind(betas, phis, etas)

# save the correponding simulation results in a dataset
Y = matrix(0, nrow = m, ncol = 4200) # the length of the data output is 4200 (100 locations, 21 time points, 2 variables)
for (s in 1:m){
  print(s)
  beta = unname(theta[s,1])
  phi = unname(theta[s,2])
  eta = unname(theta[s,3])
  
  Y[s,] = Spatial_SIR_simulation(beta,phi,eta,
                                 N,S,I,R,rows,cols,row_index, col_index, total_people,
                                 left_neighbor_I,right_neighbor_I,top_neighbor_I,bottom_neighbor_I,
                                 num_events,masking_probability,number_reporting_times,reporting_times,grid_observation_indices,
                                 dice_roll,location_dice_roll,delta_t,
                                 location_point_grid,I_grid_count,R_grid_count,
                                 aggregate_ts_plot = FALSE,lattice_snapshots = FALSE)
  

}


########################### SUMMARY STATISTICS #################################

# Use PCA to summarize the Y dataset
number_of_pcs = 3
eigen_decomp = eigen(cov(Y))
Z = Y%*%eigen_decomp$vectors[,1:number_of_pcs]

# standardize the PC scores for the neural network
mean_vector = apply(Z, 2, mean)
sd_vector = apply(Z, 2, sd)
std_Z = scale(Z)


######################### FITTING THE NEURAL NETWORK ###########################

# transform the parameters so the heterogeneous normal model makes sense for them
normal_theta = theta
normal_theta[,1] = qnorm((theta[,1] - 0.1)/.8) # beta
normal_theta[,2] = qnorm((theta[,2] - 0.1)/.8) # phi
normal_theta[,3] = qnorm((theta[,3] - 0.1)/.8) # eta

# Neural Network Architecture -- completely connected with two layers
layer1_nodes = 100
layer2_nodes = 25

# fit the neural networks (one for each parameter)
model_to_train = cbind(normal_theta, std_Z)
W  <- list()
for(j in 1:ncol(theta)){
  print(paste("Training parameter",j))
  
  gammaj     <- normal_theta[,j] # response for the jth neural network
  

  initj      <- init_hetnorm(p=ncol(std_Z),L1=layer1_nodes,L2=layer2_nodes,
                             init_mn=mean(gammaj),init_sd=sd(gammaj)) # initialize weights with normal distribution
  
  # adaptive SGD using our standardized PC Scores and the log-normal as the loss function
  modelj     <- adam(w=initj, x=std_Z, y=gammaj,
                     loss=loss_hetnorm, grad=grad_hetnorm, verbose = 1) 
  
  W[[j]] <- as.list(modelj$w) # save the weights of each neural network
}


########################### SIMULATION STUDY ###################################


# # Setting 1
# true_beta= .7
# true_phi = .8
# true_eta = .5

# Setting 2
true_beta= .5
true_phi = .3
true_eta = .3

theta0 = c(true_beta, true_phi, true_eta)

nsims = 100
q_eval = c(0.05, 0.5, 0.95) # medians for the point estimates, and 0.05, 0.95 for the 90% credible intervals

simulation_results = array(0, c(nsims, 3, length(q_eval)))
sd_results = matrix(0, nrow = nsims, ncol = 3)

for(sim in 1:nsims){
  print(sim)
  
  # generate true data from the true parameters
  Y_sampled_from_true = Spatial_SIR_simulation(true_beta,true_phi,true_eta,
                                               N,S,I,R,rows,cols,row_index, col_index, total_people,
                                               left_neighbor_I,right_neighbor_I,top_neighbor_I,bottom_neighbor_I,
                                               num_events,masking_probability,number_reporting_times,reporting_times,grid_observation_indices,
                                               dice_roll,location_dice_roll,delta_t,
                                               location_point_grid,I_grid_count,R_grid_count,
                                               aggregate_ts_plot = FALSE,lattice_snapshots = FALSE)
  
  # transfer the true data to our standardized PC scores
  true_Z = Y_sampled_from_true%*%eigen_decomp$vectors[,1:number_of_pcs]
  Z_combined = rbind(true_Z, Z)
  std_Z_combined = scale(Z_combined)
  std_true_Z = std_Z_combined[1,]
  
  # fitting the true data into our neural netowork...
  for (j in 1:length(theta0)){
    pred      <- predict_hetnorm(W[[j]],std_true_Z) # get the posterior for each parameter
    simulation_results[sim,j,] <- qnorm(q_eval,pred$mu,pred$sigma) # medians and credible intervals of each posterior
    sd_results[sim,j] = pred$sigma[1,1] # standard deviation of each posterior 
  }
}


simulation_results = pnorm(simulation_results)*.8 + .1 # transform back into uniform parameters for easier interpretation

# coverage
parameter_coverage = c(sum((simulation_results[,1,1] < true_beta) & (true_beta < simulation_results[,1,3]))/nsims,
                       sum((simulation_results[,2,1] < true_phi) & (true_phi < simulation_results[,2,3]))/nsims,
                       sum((simulation_results[,3,1] < true_eta) & (true_eta < simulation_results[,3,3]))/nsims)
print(parameter_coverage)

# MAD
parameter_MAD = c(mean(abs(simulation_results[,1,2] - true_beta)), mean(abs(simulation_results[,2,2] - true_phi)), mean(abs(simulation_results[,3,2] - true_eta)))
print(parameter_MAD)

# boxplot of posterior medians vs. true values
labs <- c(expression(beta), expression(phi), expression(eta))
c1 <- c(parameter_coverage[1], parameter_coverage[2], parameter_coverage[3])
boxplot(simulation_results[,,2],axes=F,outline=F, ylim = c(0,.9), ylab="Posterior median",cex.lab=1.25)
axis(2,cex.axis=1.25,cex.lab=1.25)
axis(1,at=1:3,lab=labs,cex.axis=1.25,cex.lab=1.25)
points(theta0,col=2,pch=19)
text(1:3,-0.05,c1,pos=3,cex=1.25)
