
# Elliot Maceda
# NCSU, May 2024

# medium-100
# 0.54 (0.54, 0.54)
# 0.96 (0.96,0.96)
# 0.14 (0.14,0.14)

# -1.18
 
###################### PACKAGES AND FUNCTIONS REQUIRED #########################
library(ggplot2)
library(MASS)
library(spam)

num_params = function(input,HL1, HL2){
  first_part = input*HL1
  second_part = HL1*HL2
  third_part = HL2*2
  result = first_part + second_part + third_part
  return(result)
}


# setwd("C:/Users/71000/Desktop/Research Projects/VaNBayes/My Work")

# Spatial SIR Simulator
source("Spatial_SIR_function_large_grids.R")

# variational posterior
source("hetnorm2.R")

# neural network optimizer
source("adam_muted.R")

# creating the grid and variables that help run the Spatial SIR simulations efficiently
source("SIR_Grid_Setup.R")

# importing this for the log-density of each population
load("zika Brazil data.rdata")
state_density = exp(log_density)
state_weights = state_density/sum(state_density)
weight_mat = matrix(state_weights, nrow = 27, ncol = 40, byrow = FALSE) # rows correspond to states
# ########## Generating NN Training Data using Training Distribution #############
# 
# set.seed(NULL)
# m <- 10000 # recommended 100000
# 
# # generate the parameters
# betas = runif(m, 0.1, .9)
# phis = runif(m, 0.1, .9)
# etas = runif(m, 0.1, .9)
# theta = cbind(betas, phis, etas)
# 
# # save the correponding simulation results in a dataset
# Y = matrix(0, nrow = m, ncol = 4200) # the length of the data output is 4200 (100 locations, 21 time points, 2 variables)
# for (s in 1:m){
#   print(s)
#   beta = unname(theta[s,1])
#   phi = unname(theta[s,2])
#   eta = unname(theta[s,3])
#   
#   Y[s,] = Spatial_SIR_simulation(beta,phi,eta,
#                                  N,S,I,R,rows,cols,row_index, col_index, total_people,
#                                  left_neighbor_I,right_neighbor_I,top_neighbor_I,bottom_neighbor_I,
#                                  num_events,masking_probability,number_reporting_times,reporting_times,grid_observation_indices,
#                                  dice_roll,location_dice_roll,delta_t,
#                                  location_point_grid,I_grid_count,R_grid_count,
#                                  aggregate_ts_plot = FALSE,lattice_snapshots = FALSE)
#   
# 
# }

########## IMPORTING THE TRAINING DATA FROM THE TRAINING DISTRIBUTION ##########

imported_data = as.matrix(read.csv("original_prior/VaNBayes_Training_Data_originalprior1.csv"))[,-1]
true_data = as.matrix(read.csv("VaNBayes_True_Data.csv"))[,-1]
Y = imported_data[,-c(1,2,3,4)]
theta = imported_data[,c(1,2,3,4)]

for (i in 2:8){
  # can copy and paste this for as many times as we have files
  more_data = as.matrix(read.csv(paste0("original_prior/VaNBayes_Training_Data_originalprior",i,".csv")))[,-1]
  more_Y = more_data[,-c(1,2,3,4)]
  more_theta = more_data[,c(1,2,3,4)]
  Y = rbind(Y, more_Y)
  theta = rbind(theta, more_theta)
}

# Y = Y[110001:120000,]
# theta = theta[110001:120000,]


dim(Y)
dim(theta)

# check if the prior is correct

true_data_mat = matrix(true_data, ncol = 40)
test = matrix(1:1080, ncol = 40)
test2 = matrix(rep(1:27, times = 40), ncol = 40)
state_populations = matrix(rep(gridpop, times = 40), ncol = 40)
state_populations_vector = c(state_populations)
true_data_mat_weighted = cbind(matrix(0, nrow = nrow(true_data_mat), ncol = 4), true_data_mat*weight_mat)
true_data_mat_counts = cbind(matrix(0, nrow = nrow(true_data_mat), ncol = 4), true_data_mat)
true_data_mat_average = apply(true_data_mat_weighted, 2, sum)
true_data_mat_total = apply(true_data_mat_counts, 2, sum)


true_data_state_totals = apply(true_data_mat, 1, sum)
true_data_state_time_of_totals = apply(true_data_mat, 1, which.max)


total_reps = nrow(Y)
most_populous_mat = matrix(0, nrow = total_reps, ncol = 40+4) # 40 times, 4 parameters to make
most_dense_mat = matrix(0, nrow = total_reps, ncol = 40+4)
average_mat = matrix(0, nrow = total_reps, ncol = 40+4)
total_mat = matrix(0, nrow = total_reps, ncol = 40+4)


state_totals_mat = matrix(0, nrow = total_reps, ncol = 27)
state_time_of_totals_mat = matrix(0, nrow = total_reps, ncol = 27)

for (i in 1:total_reps){
  Y_mat = matrix(Y[i,], ncol = 40)
  Y_mat_average = apply(Y_mat*weight_mat, 2, sum)
  Y_mat_total = apply(Y_mat, 2, sum)
  most_populous_mat[i,] = c(theta[i,],Y_mat[25,])
  most_dense_mat[i,] = c(theta[i,],Y_mat[7,])
  average_mat[i,] = c(theta[i,],Y_mat_average)
  total_mat[i,] = c(theta[i,],Y_mat_total)
  
  # used for summary statistics (additional covariates to help)
  state_totals_mat[i,] = apply(Y_mat, 1, sum)
  state_time_of_totals_mat[i,] = apply(Y_mat, 1, which.max)
}

# most_populous_matplot = rbind(true_data_mat[25,], most_populous_mat)
# most_dense_matplot = rbind(true_data_mat[7,], most_dense_mat)
average_matplot = rbind(true_data_mat_average, average_mat)
total_matplot = rbind(true_data_mat_total, total_mat)


# # Euclidean Distance Calculation
# # sweep() - performs an operation across a dimension of a matrix
# subtracted_populous = sweep(most_populous_matplot[-1, -c(1:4)], 2, most_populous_matplot[1, -c(1:4)], "-") #column-wise, subtract our real data
# Euclidean_distance_populous = apply(subtracted_populous^2, 1, sum) # entrywise we square, then rowwise we add
# most_similar_populous_params = most_populous_matplot[which.min(Euclidean_distance_populous) + 1,c(1:4)] # row of the closest observation
# most_similar_populous_params
# 
# subtracted_dense = sweep(most_dense_matplot[-1, -c(1:4)], 2, most_dense_matplot[1, -c(1:4)], "-") #column-wise, subtract our real data
# Euclidean_distance_dense = apply(subtracted_dense^2, 1, sum) # entrywise we square, then rowwise we add
# most_similar_dense_params = most_dense_matplot[which.min(Euclidean_distance_dense) + 1,c(1:4)] # row of the closest observation
# most_similar_dense_params
# 
# subtracted_average = sweep(average_matplot[-1, -c(1:4)], 2, average_matplot[1, -c(1:4)], "-") #column-wise, subtract our real data
# Euclidean_distance_average = apply(subtracted_average^2, 1, sum) # entrywise we square, then rowwise we add
# most_similar_average_params = average_matplot[which.min(Euclidean_distance_average) + 1,c(1:4)] # row of the closest observation
# most_similar_average_params



matplot(t(true_data_mat[,-c(1:4)]), type = "l", lty = 5, lwd = 3)
title("true data")
matplot(t(average_matplot[,-c(1:4)]),
        type = "l",
        col = c("red", rep("gray", total_reps)),
        lty = c(1, rep(3, total_reps)),
        lwd = c(3, rep(1, total_reps)),
        xlab = "Week",
        ylab = "Average Counts of Infected, Weighted by State Density") # average counts or total countS?
# title("average cases per state")


integrated_avg_trajectories = apply(average_matplot[-1,-c(1:4)], 1, sum)
true_avg_integral = sum(average_matplot[1,-c(1:4)])
hist(integrated_avg_trajectories, xlab = "Integral of Trajectories", main = "")
abline(v = true_avg_integral, col = "red", lwd = 3)
# Y = sweep(Y, 2, FUN = "/", state_populations_vector)
eigen_decomp = eigen(cov(Y))

# integrated_total_trajectories = apply(total_matplot[-1,-c(1:4)], 1, sum)
# true_total_integral = sum(total_matplot[1,-c(1:4)])
# hist(integrated_total_trajectories, xlab = "Integral of Trajectories", main = "")
# abline(v = true_total_integral, col = "red", lwd = 3)

########################### SUMMARY STATISTICS #################################

# Neural Network Architecture -- completely connected with two layers
neural_network_size = "medium"
summary_vars = 20 # 4 (83%), 6 (90%), 9 (95%), 20 (99%), 30 (99.5%)
# PC 5
# 15
# 40 
# 94 

# beta0_lower = -0.9
# beta0_upper = 0.1
# beta1_lower = 0
# beta1_upper = 0.7
# phi_lower = 0
# phi_upper = 0.24
# nu_lower = 1.01
# nu_upper = 10

beta0_lower = -3
beta0_upper = 1
beta1_lower = -1
beta1_upper = 1
# phi_lower = -2
# phi_upper = 1
nu_lower = 1.01
nu_upper = 10


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

# if (summary_vars == 5){
#   # # Use PCA to summarize the Y dataset
#   number_of_pcs = 5
#   eigen_decomp = eigen(cov(Y))
#   Z = Y%*%eigen_decomp$vectors[,1:number_of_pcs]
#   true_Z = true_data%*%eigen_decomp$vectors[,1:number_of_pcs]
# } else if (summary_vars == 15){
#   # # Use PCA to summarize the Y dataset
#   number_of_pcs = 15
#   eigen_decomp = eigen(cov(Y))
#   Z = Y%*%eigen_decomp$vectors[,1:number_of_pcs]
#   true_Z = true_data%*%eigen_decomp$vectors[,1:number_of_pcs]
# } else if (summary_vars == 20){
#   # # Use PCA to summarize the Y dataset
#   number_of_pcs = 20
#   eigen_decomp = eigen(cov(Y))
#   Z = Y%*%eigen_decomp$vectors[,1:number_of_pcs]
#   true_Z = true_data%*%eigen_decomp$vectors[,1:number_of_pcs]
# } else if (summary_vars == 50){
#   # # Use PCA to summarize the Y dataset
#   number_of_pcs = 50
#   eigen_decomp = eigen(cov(Y))
#   Z = Y%*%eigen_decomp$vectors[,1:number_of_pcs]
#   true_Z = true_data%*%eigen_decomp$vectors[,1:number_of_pcs]
# } else if (summary_vars == 100){
#   # # Use PCA to summarize the Y dataset
#   number_of_pcs = 100
#   eigen_decomp = eigen(cov(Y))
#   Z = Y%*%eigen_decomp$vectors[,1:number_of_pcs]
#   true_Z = true_data%*%eigen_decomp$vectors[,1:number_of_pcs]
# } else if (summary_vars == 200){
#   # # Use PCA to summarize the Y dataset
#   number_of_pcs = 200
#   eigen_decomp = eigen(cov(Y))
#   Z = Y%*%eigen_decomp$vectors[,1:number_of_pcs]
#   true_Z = true_data%*%eigen_decomp$vectors[,1:number_of_pcs]
# } else if (summary_vars == 40){
#   # Averaged counts wrt density
#   Z = average_matplot[-1,-c(1:4)]
#   true_Z = average_matplot[1,-c(1:4)]
# } else {
#   # Averaged counts wrt density
#   Z = average_matplot[-1,-c(1:4)]
#   true_Z = average_matplot[1,-c(1:4)]
#   # additional predictors to help:
#   Z = cbind(Z,state_totals_mat,state_time_of_totals_mat)
#   true_Z = c(true_Z, true_data_state_totals,true_data_state_time_of_totals)
# }


# # Use PCA to summarize the Y dataset
number_of_pcs = summary_vars
# eigen_decomp = eigen(cov(Y))
Z = Y%*%eigen_decomp$vectors[,1:number_of_pcs]
true_Z = true_data%*%eigen_decomp$vectors[,1:number_of_pcs]

reconstructed_true_data = true_Z%*%t(eigen_decomp$vectors[,1:number_of_pcs])
reconstructed_true_data_mat = matrix(reconstructed_true_data, ncol = 40)
matplot(t(reconstructed_true_data_mat), type = "l", lty = 5, lwd = 3)
title("reconstructed true data")
matplot(t(true_data_mat), type = "l", lty = 5, lwd = 3)
title("true data")


num_PCs = 90

reconstructed_true_data = true_data%*%eigen_decomp$vectors[,1:num_PCs]%*%t(eigen_decomp$vectors[,1:num_PCs])
reconstructed_true_data_mat = matrix(reconstructed_true_data, ncol = 40)
matplot(t(true_data_mat), type = "l", lty = 5, lwd = 3, ylab = "Counts of Infected", xlab = "Weeks")
title("true data")
matplot(t(reconstructed_true_data_mat), type = "l", lty = 5, lwd = 3, ylab = "Counts of Infected", xlab = "Weeks")
title(paste0("reconstructed true data with ", num_PCs, " PC's"))


# standardize the PC scores for the neural network
Z_combined = rbind(true_Z, Z)
std_Z_combined = scale(Z_combined)
std_true_Z = std_Z_combined[1,]
std_Z = std_Z_combined[-1,]


######################### FITTING THE NEURAL NETWORK ###########################


# These are going to be pi(theta)/Pi(theta)
myprior_weights = matrix(0, nrow = nrow(theta), ncol = ncol(theta))

# These are Parker's hyperparameters (the model's prior)
myprior_weights[,1] = dunif(theta[,1], -3, 1)# pi(theta)
myprior_weights[,2] = dunif(theta[,2], -1, 1)# pi(theta)
myprior_weights[,3] = dlnorm(theta[,3], -2, 1)# pi(theta)
myprior_weights[,4] = dunif(theta[,4], 1.01, 10)# pi(theta)

for(j in 1:ncol(theta)){ # Dividing by Pi(theta)
  # how to get Pi for each parameter
  d <- density(theta[,j])
  dd <- approxfun(d$x, d$y)
  trainingDensity = dd(theta[,j])
  myprior_weights[,j] <- myprior_weights[,j]/trainingDensity
}





# transform the parameters so the heterogeneous normal model makes sense for them
normal_theta = theta
normal_theta[,1] = qnorm((theta[,1] - (beta0_lower))/(beta0_upper-(beta0_lower)))       # beta0
normal_theta[,2] = qnorm((theta[,2] - (beta1_lower))/(beta1_upper-(beta1_lower)))       # beta1
# normal_theta[,3] = qnorm((theta[,3] - (phi_lower))/(phi_upper-(phi_lower)))       # phi
normal_theta[,3] = log(theta[,3])       # phi
normal_theta[,4] = qnorm((theta[,4] - (nu_lower))/(nu_upper-(nu_lower)))


# fit the neural networks (one for each parameter)



W  <- list()
for(j in 1:ncol(theta)){
  print(paste("Training parameter",j))
  
  gammaj     <- normal_theta[,j]# theta[,j]# normal_theta[,j] # response for the jth neural network
  

  initj      <- init_hetnorm(p=ncol(std_Z),L1=layer1_nodes,L2=layer2_nodes,
                             init_mn=mean(gammaj),init_sd=sd(gammaj)) # initialize weights with normal distribution
  
  # adaptive SGD using our standardized PC Scores and the log-normal as the loss function
  modelj     <- adam(w=initj, x=std_Z, y=gammaj, priorweights = myprior_weights[,j], lr = 1e-2, early_stop = 200,
                     loss=loss_hetnorm, grad=grad_hetnorm, verbose = 1, batchsize = 1000) 
  
  W[[j]] <- as.list(modelj$w) # save the weights of each neural network
}


########################### Estimating the Posterior ###########################



# fitting the true data into our neural netowork...
pred      <- predict_hetnorm(W[[1]],std_true_Z)
normal_posterior = rnorm(10000,mean = pred$mu, sd = pred$sigma)
# plot(density(normal_posterior))
beta0_posterior = pnorm(normal_posterior)*(beta0_upper-(beta0_lower)) + (beta0_lower)
plot(density(beta0_posterior))
lines(density(theta[,1]), col = "red")


pred      <- predict_hetnorm(W[[2]],std_true_Z)
pred$mu
normal_posterior = rnorm(10000,mean = pred$mu, sd = pred$sigma)
# plot(density(normal_posterior))
beta1_posterior = pnorm(normal_posterior)*(beta1_upper-(beta1_lower)) + (beta1_lower)
plot(density(beta1_posterior))
lines(density(theta[,2]), col = "red")

pred      <- predict_hetnorm(W[[3]],std_true_Z)
pred$mu
normal_posterior = rnorm(10000,mean = pred$mu, sd = pred$sigma)
# plot(density(normal_posterior))
# phi_posterior = pnorm(normal_posterior)*(phi_upper-(phi_lower)) + (phi_lower)
phi_posterior = exp(normal_posterior)
plot(density(phi_posterior))
lines(density(theta[,3]), col = "red")

pred      <- predict_hetnorm(W[[4]],std_true_Z)
pred$mu
normal_posterior = rnorm(10000,mean = pred$mu, sd = pred$sigma)
# plot(density(normal_posterior))
nu_posterior = pnorm(normal_posterior)*(nu_upper-(nu_lower)) + (nu_lower)
plot(density(nu_posterior))
lines(density(theta[,4]), col = "red")
# print(c(mean(nu_posterior), quantile(nu_posterior, 0.025), quantile(nu_posterior, 0.975)))

# Fits
# NOTE: we use std_Z and normal_theta here (originating from theta)
# however, in reality these need to be generated independently
LS  <- rep(0,4)
PIT <- matrix(0,nrow(std_Z),4)
for(j in 1:4){
  pred    <- predict_hetnorm(W[[j]],std_Z)
  LS[j]   <- mean(dnorm(normal_theta[,j],pred$mu,pred$sigma,log=TRUE))
  PIT[,j] <- pnorm(normal_theta[,j],pred$mu,pred$sigma)
}
png(paste0("final_results2/", neural_network_size, "_", summary_vars, "_PIT.png"), width = 700, height = 700, pointsize = 20) # can also use "bmp", "jpeg", "tiff", etc.
plot(NA,xlim=0:1,ylim=0:1,xlab="Theoretical quantile",ylab="Observed quantile")
mycolors = c("red", "blue", "green")
for(j in 1:3){
  qqq   <- seq(0.01,0.99,0.01)
  lines(qqq,quantile(PIT[,j],qqq),lty=j, col = mycolors[j])
}
legend(x = c(0.8,1.00), y= c(0,0.2), legend = c(expression(beta[0]), expression(beta[1]), expression(phi)), col = mycolors, lty = c(1,2,3))
dev.off()
print(c(mean(beta0_posterior), quantile(beta0_posterior, 0.025), quantile(beta0_posterior, 0.975)))
print(c(mean(beta1_posterior), quantile(beta1_posterior, 0.025), quantile(beta1_posterior, 0.975)))
print(c(mean(phi_posterior), quantile(phi_posterior, 0.025), quantile(phi_posterior, 0.975)))
LS
mean(LS)
print(neural_network_size)
print(summary_vars)




############################ IMPORTANCE ESTIMATES ##############################

library(ALEPlot)

predict_for_ALEPlot_mu = function(X.model, newdata){ # newdata is 130000 by 40... want 130000 different predictions
  
  predict_hetnorm_numeric = function(X, mymodel){
    return(predict_hetnorm(mymodel, X)$mu[1,1])
  }
  
  results = apply(X = newdata, MARGIN = 1, FUN = predict_hetnorm_numeric, X.model)
  return(results)
}

predict_for_ALEPlot_sigma = function(X.model, newdata){ # newdata is 130000 by 40... want 130000 different predictions
  
  predict_hetnorm_numeric = function(X, mymodel){
    return(predict_hetnorm(mymodel, X)$sigma[1,1])
  }
  
  results = apply(X = newdata, MARGIN = 1, FUN = predict_hetnorm_numeric, X.model)
  return(results)
}

j = 1
ALEPlot(X = std_Z, X.model = W[[1]], pred.fun = predict_for_ALEPlot_mu, J = j) # horizontal is no importance... movement is importance
j = 2
par(new = TRUE)
ALEPlot(X = std_Z, X.model = W[[1]], pred.fun = predict_for_ALEPlot_mu, J = j) # horizontal is no importance... movement is importance





########################## ALE PLOTS FOR MU ####################################

# make these like the PIT plots--- different colors for each variable


ALE_mycolors = c("red", "blue", "green")


num_covs = 3
num_top = 2

ALE_importance_mat = matrix(0, nrow = num_covs, ncol = 3)
rownames(ALE_importance_mat) = as.character(1:num_covs) # can change this to be more specific
saved_ALE_plots = vector(mode = "list", length = 3)
names(saved_ALE_plots) = c("beta0","beta1","phi")

for(j in 1:3){ # beta0, beta1, phi
  saved_ALE_plots[[j]] = vector(mode = "list", length = num_covs)
  for(i in 1:num_covs){ # covariates
    ALE_obj = ALEPlot(X = std_Z, X.model = W[[j]], pred.fun = predict_for_ALEPlot_mu, J = i)
    saved_ALE_plots[[j]][[i]] = ALE_obj
    importance_metric = (ALE_obj$f.values[which.max(ALE_obj$x.values)] - ALE_obj$f.values[which.min(ALE_obj$x.values)])/(max(ALE_obj$x.values) - min(ALE_obj$x.values))
    ALE_importance_mat[i,j] = importance_metric
  }
}


for(j in 1:3){ # beta0, beta1, phi
  topn_cov_indices = sort(abs(ALE_importance_mat[,j]), index.return = TRUE)$ix[1:num_top]
  plot(NA,xlim=c(-10,10),ylim=c(-0.5,0.5),xlab="x-value",ylab="ALE_plot value", main = names(saved_ALE_plots)[j])
  for(i in 1:num_top){ # top n covariates
    ALE_plot_info = saved_ALE_plots[[j]][[topn_cov_indices[i]]]
    lines(ALE_plot_info$x.values, ALE_plot_info$f.values, col = ALE_mycolors[i], lty = 1:num_top)
  }
  legend(x = "bottomright", legend = rownames(ALE_importance_mat)[topn_cov_indices], col = ALE_mycolors[1:num_top], lty = 1:num_top)
}


########################## ALE PLOTS FOR SIGMA #################################

ALE_importance_mat = matrix(0, nrow = 40, ncol = 3)
saved_ALE_plots = vector(mode = "list", length = 3)
names(saved_ALE_plots) = c("beta0","beta1","phi")

for(j in 1:3){ # beta0, beta1, phi
  saved_ALE_plots[[j]] = vector(mode = "list", length = 40)
  for(i in 1:40){ # covariates
    ALE_obj = ALEPlot(X = std_Z, X.model = W[[j]], pred.fun = predict_for_ALEPlot_mu, J = i)
    saved_ALE_plots[[j]][[i]] = ALE_obj
    importance_metric = (ALE_obj$f.values[which.max(ALE_obj$x.values)] - ALE_obj$f.values[which.min(ALE_obj$x.values)])/(max(ALE_obj$x.values) - min(ALE_obj$x.values))
    ALE_importance_mat[i,j] = importance_metric
  }
}


for(j in 1:3){ # beta0, beta1, phi
  top5cov_indices = sort(abs(ALE_importance_mat[,j]), index.return = TRUE)$ix[1:5]
  for(i in 1:length(top5cov_indices)){ # top 5 covariates
    ALE_plot_info = saved_ALE_plots[[j]][[top5cov_indices[i]]]
    lines(ALE_plot_info$x.values, ALE_plot_info$f.values, col = mycolors[i])
  }
}



############################ PCA IMAGE PICS ####################################

# take the first row of the loading matrix (corresponds to a particular time/space)
# loading = eigen_decomp$vectors[,1:number_of_pcs]
# compare image(matrix(loading[,1], nrow = 40, ncol = 27, byrow = FALSE)) to image(true_data_mat)
# order the spatial columns/rows by population density or latitude/longitude

# true_data # this is just imported
# true_data_mat = matrix(true_data, ncol = 40, byrow = F) # this is how we transform it to a matrix we are happy with
true_data_mat_ordered = true_data_mat[order(log_density, decreasing = T), ]

true_data_mat_rotated <- apply(true_data_mat_ordered, 2, rev)
image(1:ncol(true_data_mat_ordered), 1:nrow(true_data_mat_ordered), t(true_data_mat_rotated), ylab = "State (from most to least dense)", xlab = "Week", yaxt = "n")

loading_matrix = eigen_decomp$vectors[,1:number_of_pcs]

for (k in 1:number_of_pcs){
  loading_weights = loading_matrix[,k]
  loading_weights_mat = matrix(loading_weights, ncol = 40, byrow = F)
  loading_weights_mat_ordered = loading_weights_mat[order(log_density, decreasing = T), ]
  
  loading_weights_mat_ordered
  true_data_mat_rotated <- apply(loading_weights_mat_ordered, 2, rev)
  image(1:ncol(loading_weights_mat_ordered), 1:nrow(loading_weights_mat_ordered), t(true_data_mat_rotated), ylab = "State", xlab = "Week", yaxt = "n")
}

# find the most important variables for each parameter/variational parameter and interpret

order(log_density, decreasing = T)

boxplot(loading_matrix)

# the 2nd most dense is 19 (Rio de Janeiro, 3rd most populous, sourthern part... behind Sao Paulo but close to it geographically)
# the 15th most dense is 5 (Bahia, 4th most populous, in the more northern part... touches 9 states!)
# the 25th most dense is 11 (Mato Grosso, 18th most populous, western part of Brazil, touches 6 states... maybe "connects" the states in the amazon rain forest to the rest of Brazil (at least one of the states that does this))
# honorable mentions: 
# Distrito Federal (most dense, but located within Goias and tiny population)
# Sao Paulo: third most dense and 1st most populous, but close to Rio... so I guess the PCA prefers Rio because of its greater density
# Amazonas: Also might be pretty important to PCA? The NW-most state with a fair amount of people and borders 5 states (second least dense)
# Minas Gerais: 1-more dense than Bahia, located between Bahia and Rio/Sao Paulo



############################### MAP OF BRAZIL ##################################


# install.packages("geobr")
library(geobr)
library(ggplot2)
# Read all municipalities in the country at a given year
states = read_state(code_state="all", year=2016)

states = states %>% 
  arrange(name_state)

pop_density = exp(log_density)
names(pop_density) = "Population Density"
# Map of Brazil for PCA analysis
fill_colors_PCA = case_when(states$name_state == "Rio De Janeiro" ~ "red",
                        states$name_state == "Bahia" ~ "orange",
                        states$name_state == "Mato Grosso" ~ "blue",
                        TRUE ~ "white")

png(paste0("results/Brazil_Map_PCA.png"), width = 700, height = 700, pointsize = 20) # can also use "bmp", "jpeg", "tiff", etc.
ggplot(states) + geom_sf(fill = fill_colors_PCA, color = "black") + theme_void()
dev.off()

png(paste0("results/Brazil_Map_Density.png"), width = 700, height = 700, pointsize = 20) # can also use "bmp", "jpeg", "tiff", etc.
ggplot(states) + geom_sf(aes(fill = pop_density), color = "black") + theme_void() + guides(fill = guide_legend(title="Population Density")) 
dev.off()

ggplot(states) + geom_sf(aes(fill = pop_density), color = "black") + theme_void() + labs(fill = "Population Density") 
ggplot(states) + geom_sf(aes(fill = pop_density), color = "black") + theme_void() + scale_fill_continuous("Population Density") 

################################### SUMMARY STATISTICS INTERPRETATION ##########
# Summary Statistic list: 


sum.stat.names1 = paste0("Total in Week ", 1:40)
sum.stat.names2 = paste0("Total in ", states_names)
sum.stat.names3 = paste0("When Infections Peak in ", states_names)

sum.stat.names = c(sum.stat.names1, sum.stat.names2, sum.stat.names3)

# MU

# beta0 -- all fairly strong negative effects
sum.stat.names[2]
sum.stat.names[32]
sum.stat.names[3]
sum.stat.names[1]
sum.stat.names[48]

# beta1 -- all fairly strong positive effects
sum.stat.names[2]
sum.stat.names[1]
sum.stat.names[3]
sum.stat.names[10]
sum.stat.names[15]

# phi -- 2 is positive, the rest are negative effects
sum.stat.names[61]
sum.stat.names[54]
sum.stat.names[2]
sum.stat.names[65]
sum.stat.names[50]

fill_colors_MU = case_when(states$name_state == "São Paulo" ~ "red",
                           states$name_state == "Rio Grande Do Sul" ~ "orange",
                           states$name_state == "Espírito Santo" ~ "blue",
                           states$name_state == "Maranhão" ~ "purple",
                              TRUE ~ "white")
png(paste0("results/Brazil_Map_MU.png"), width = 700, height = 700, pointsize = 20) # can also use "bmp", "jpeg", "tiff", etc.
ggplot(states) + geom_sf(fill = fill_colors_MU, color = "black") + theme_void()
dev.off()

# SIGMA

# beta0 -- 4,3,2 are much more important than the other times
sum.stat.names[4]
sum.stat.names[3]
sum.stat.names[2]
sum.stat.names[15]
sum.stat.names[8]

# beta1 -- all important and proportional
sum.stat.names[3]
sum.stat.names[40]
sum.stat.names[32]
sum.stat.names[4]
sum.stat.names[2]

# phi
sum.stat.names[65] # only this one... the others aren't important (directly proportional)
sum.stat.names[88]
sum.stat.names[1]
sum.stat.names[40]
sum.stat.names[61]


fill_colors_SIGMA = case_when(states$name_state == "São Paulo" ~ "red",
                           states$name_state == "Rio Grande Do Sul" ~ "orange",
                           TRUE ~ "white")

png(paste0("results/Brazil_Map_SIGMA.png"), width = 700, height = 700, pointsize = 20) # can also use "bmp", "jpeg", "tiff", etc.
ggplot(states) + geom_sf(fill = fill_colors_SIGMA, color = "black") + theme_void()
dev.off()

######################### FACETED ERROR BARS (DEMO) ###################################


library(ggplot2)
df <- data.frame(
  trt = factor(c(1,2,3,4)),
  resp = c(1,5,3,4),
  upper = c(1.1,5.3, 3.3, 4.2),
  lower = c(0.8, 4.6, 2.4, 3.6)
)

factor1 = factor(c("level 1", "level 1", "level 1", "level 1"))
factor2 = factor(c("A", "A", "A", "A")) 
df1 = cbind(df, factor1, factor2)
factor1 = factor(c("level 2", "level 2", "level 2", "level 2"))
factor2 = factor(c("A", "A", "A", "A")) 
df2 = cbind(df, factor1, factor2)
factor1 = factor(c("level 1", "level 1", "level 1", "level 1"))
factor2 = factor(c("B", "B", "B", "B")) 
df3 = cbind(df, factor1, factor2)
factor1 = factor(c("level 2", "level 2", "level 2", "level 2"))
factor2 = factor(c("B", "B", "B", "B")) 
df4 = cbind(df, factor1, factor2)

big_df = rbind(df1, df2, df3, df4)
levels(big_df$factor1) = c(expression(alpha), expression(beta))
# levels(big_df$factor1) = c(expression("100 "*.L^"-1"*""), expression("200 "*.L^"-1"*""))

p <- ggplot(big_df, aes(trt, resp)) 

p + geom_point(aes(x = trt, y = resp)) + geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) + facet_grid(factor1~factor2, labeller = label_parsed) + ylab("") + xlab("x axis label")

######################### FACETED ERROR BARS ###################################


library(ggplot2)


errorbars_df = expand.grid(PCs = factor(c(4,6,9,20)),
                           parameter = factor(c("beta0","beta1","phi")),
                           NN_size = factor(c("Small", "Medium", "Large"), levels = c("Small", "Medium", "Large")))
posterior_medians = c(-0.32, -0.47, -0.62, -1.16, 0.27, 0.56, 0.44, 0.51, 0.17, 0.12, 0.10, 0.15, -0.46, -0.67, -0.53, -1.33, 0.19, 0.27, 0.38, 0.60, 0.09, 0.10, 0.14, 0.13, -0.25, -0.53, -0.55, -0.78, 0.17, 0.26, 0.35, 0.71, 0.18, 0.08, 0.09, 0.23)
upper_credible = c(-1.01, -1.85, -1.04, -1.18, -0.17, 0.27, 0.40, -0.54, 0.02, 0.03, 0.10, 0.13, -0.98, -1.80, -0.70, -1.47, -0.15, -0.10, 0.35, 0.14, 0.01, 0.01, 0.13, 0.11, -0.90, -1.65, -0.57, -0.80, -0.29, -0.16, 0.30, 0.61, 0.00, 0.01, 0.09, 0.18)
lower_credible = c(0.27, 0.64, -0.22, -1.13, 0.64, 0.78, 0.48, 0.99, 0.71, 0.30, 0.10, 0.17, 0.04, 0.35, -0.37, -1.19, 0.51, 0.62, 0.41, 0.89, 0.34, 0.40, 0.15, 0.17, 0.31, 0.42, -0.53, -0.75, 0.60, 0.63, 0.41, 0.80, 1.05, 0.31, 0.10, 0.29)

errorbars_df = data.frame(errorbars_df, posterior_medians, upper_credible, lower_credible)

levels(errorbars_df$parameter) = c(expression(beta[0]), expression(beta[1]),expression(phi))
levels(errorbars_df$NN_size) = c(expression(paste("(50", ",", "10)")), expression(paste("(100", ",", "25)")),expression(paste("(150", ",", "50)")))

p <- ggplot(errorbars_df, aes(PCs, posterior_medians)) 

p + geom_point(aes(x = PCs, y = posterior_medians)) + geom_hline(yintercept = 0, color = "darkgray") + geom_errorbar(aes(ymin = lower_credible, ymax = upper_credible), width = 0.2) + facet_grid(parameter~NN_size, labeller = label_parsed, scales = "free_y") + ylab("") + xlab("Number of PC's")

# expression("100 µg " * .L^"-1" * "")
