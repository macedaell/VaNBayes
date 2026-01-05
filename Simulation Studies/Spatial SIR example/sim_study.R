
# Elliot Maceda
# NCSU, May 2024
# Spatial SIR Simulation Study

###################### PACKAGES AND FUNCTIONS REQUIRED #########################
library(ggplot2)
library(MASS)
library(spam)
library(doParallel)

# setwd("C:/")

# Spatial SIR Simulator
source("Spatial_SIR_function_large_grids.R")

# neural network optimizer
source("./adam.R")

# creating the grid and variables that help run the Spatial SIR simulations efficiently
I0 = 10
set.seed(1) # setting this to be constant so we get the same grid each time (later we'll call set.seed(NULL) to undo this)
source("SIR_Grid_Setup.R")

########## Generating NN Training Data using Training Distribution #############

# set.seed(NULL)
m <- 10000 # recommended 100000
number_of_pcs = 950 # 3 (50%), 370 (70%), 950 (90%)

# generate the parameters
peak_counts = c()
I_t = matrix(0, nrow = m, ncol = num_events)
# total_counts = c()
# peak_time = c()
# eradication_time = c()
# halftime_infection_count = c()
# save the correponding simulation results in a dataset
Y = matrix(0, nrow = m, ncol = 4200) # the length of the data output is 4200 (100 locations, 21 time points, 2 variables)
for (s in 1:m){
  print(s)
  set.seed(s + 1)
  beta = runif(1, 0.1, .9)
  phi = runif(1, 0.1, .9)
  eta = runif(1, 0.1, .9)

  # source("SIR_Grid_Setup.R")

  sim_output = Spatial_SIR_simulation(beta,phi,eta,
                                 N,S,I,R,rows,cols,row_index, col_index, total_people,
                                 left_neighbor_I,right_neighbor_I,top_neighbor_I,bottom_neighbor_I,
                                 num_events,masking_probability,number_reporting_times,reporting_times,grid_observation_indices,
                                 dice_roll,location_dice_roll,delta_t,
                                 location_point_grid,I_grid_count,R_grid_count,
                                 aggregate_ts_plot = FALSE,lattice_snapshots = FALSE)
  
  Y[s,] = sim_output$response
  
  peak_counts = c(peak_counts,sim_output$I_summary)
  I_t[s,] = sim_output$I_trajectory
  # total_counts = c(total_counts,sim_output$I_summary[2])
  # peak_time = c(peak_time,sim_output$I_summary[3])
  # eradication_time = c(eradication_time,sim_output$I_summary[4])
  # halftime_infection_count = c(halftime_infection_count, sim_output$I_summary[5])
}
I_summaries = list(peak_counts,I_t)

# hist(I_summaries[[1]]) # probably these are the best 
# hist(I_summaries[[2]]) # 
# hist(I_summaries[[3]])
# hist(I_summaries[[4]])
# hist(I_summaries[[5]])

########################### SUMMARY STATISTICS #################################

# Use PCA to summarize the Y dataset
eigen_decomp = eigen(cov(Y))
Z = Y%*%eigen_decomp$vectors[,1:number_of_pcs]

# standardize the PC scores for the neural network
mean_vector = apply(Z, 2, mean)
sd_vector = apply(Z, 2, sd)
std_Z = scale(Z)

# save the results so we don't need to run this again
save(Y, file = "Y.RData")
save(eigen_decomp, file = "eigen_decomp.RData")
save(std_Z, file = "std_Z.RData")
save(I_summaries, file = "I_summaries.RData")

######################### FITTING THE NEURAL NETWORK ###########################

load("Y.RData")
load("eigen_decomp.RData")
load("std_Z.RData")
load("I_summaries.RData")

# VaNBayes settings

# Things that probably won't change:
EPOCHS = 200
PROPVAL = 0.2
EARLY_STOP = 20

# Data settings: 
RESPONSE = .3 # choices: 1 or num between 0 and 1

# Neural Network settings:
NB_versions = c(2,2)
ARCHITECTUREs = list(c(100),c(70))
LRs = c(0.0005,0.0005) # 0.0005
BATCHSIZEs = c(100,100)  #


if (RESPONSE == 1) gammaj = I_summaries[[1]] else gammaj = I_summaries[[2]][,round(quantile(1:num_events,RESPONSE))]
if (RESPONSE == 1) choice_of_I = 1 else choice_of_I = 2
if (NB_versions[choice_of_I] == 1){
  
  source("./NegativeBinomial1.R")
  
  MOM_r = mean(gammaj)/(var(gammaj) - mean(gammaj))
  MOM_p = mean(gammaj)/var(gammaj)
  
  # NB1
  init <- init_NegativeBinomial1(p = ncol(std_Z),
                                 hidden_layers = ARCHITECTUREs[[choice_of_I]],
                                 init_r = MOM_r,
                                 init_p = MOM_p)
  
  model <- adam(w = init,
                x = std_Z,
                y = gammaj, 
                loss = loss_NegativeBinomial1,
                grad = grad_NegativeBinomial1,
                batchsize = BATCHSIZEs[choice_of_I],
                epochs = EPOCHS,
                propval = PROPVAL, 
                early_stop = EARLY_STOP,
                lr = LRs[choice_of_I],
                verbose = 2)
  W <- as.list(model$w)
  
  
  predict_function = predict_NegativeBinomial1
  
  obtain_quantiles_and_mean = function(pred, q_eval) c(qnbinom(q_eval, size = pred$r, prob = pred$p), pred$r*(1-pred$p)/pred$p)
  obtain_rank_prob = function(pred, parameter) pnbinom(parameter, size = pred$r, prob = pred$p)
  
} else { # version 2
  
  source("./NegativeBinomial2.R")
  
  init <- init_NegativeBinomial2(p = ncol(std_Z),
                                 hidden_layers = ARCHITECTUREs[[choice_of_I]],
                                 init_mn = mean(gammaj))
  
  model <- adam(w = init,
                x = std_Z,
                y = gammaj, 
                loss = loss_NegativeBinomial2,
                grad = grad_NegativeBinomial2,
                batchsize = BATCHSIZEs[choice_of_I],
                epochs = EPOCHS,
                propval = PROPVAL, 
                early_stop = EARLY_STOP,
                lr = LRs[choice_of_I],
                verbose = 1)
  W <- as.list(model$w)
  
  
  predict_function = predict_NegativeBinomial2
  
  obtain_quantiles_and_mean = function(pred, q_eval) c(qnbinom(q_eval, mu = pred$mu, size = pred$psi), pred$mu)
  obtain_rank_prob = function(pred, parameter) pnbinom(parameter, mu = pred$mu, size = pred$psi)
}

########################### SIM STUDY ###################################


# # Setting 1
# true_beta= .7
# true_phi = .8
# true_eta = .5

m_val = 10000

# Setting 2
true_beta= .5
true_phi = .3
true_eta = .3

# save the generated true I summaries 
true_peak_counts = c()
true_total_counts = c()
true_peak_time = c()
true_eradication_time = c()

# qeval
q_eval = c(0.05, 0.5, 0.95)

# simulation_results = array(0, c(m, length(q_eval)+1))
# unif_rv = numeric(m)
# for (s in 1:m){
cl = makeCluster(detectCores()-1)
registerDoParallel(cl)
simulation_results = foreach(s = 1:m_val, .combine=rbind)%dopar%{
  print(s)
  set.seed(s + m + 1)
  # source("SIR_Grid_Setup.R")
  
  sim_output = Spatial_SIR_simulation(true_beta,true_phi,true_eta,
                                      N,S,I,R,rows,cols,row_index, col_index, total_people,
                                      left_neighbor_I,right_neighbor_I,top_neighbor_I,bottom_neighbor_I,
                                      num_events,masking_probability,number_reporting_times,reporting_times,grid_observation_indices,
                                      dice_roll,location_dice_roll,delta_t,
                                      location_point_grid,I_grid_count,R_grid_count,
                                      aggregate_ts_plot = FALSE,lattice_snapshots = FALSE)

  true_peak_counts = c(true_peak_counts,sim_output$I_summary)
  
  # true_peak_counts = c(true_peak_counts,sim_output$I_summary[[1]])
  # true_total_counts = c(true_total_counts,sim_output$I_summary[2])
  # true_peak_time = c(true_peak_time,sim_output$I_summary[3])
  # true_eradication_time = c(true_eradication_time,sim_output$I_summary[4])
  
  Y_sampled_from_true = sim_output$response

  # transfer the true data to our standardized PC scores
  true_Z = Y_sampled_from_true%*%eigen_decomp$vectors[,1:number_of_pcs]
  
  std_true_Z = (true_Z - attributes(std_Z)$`scaled:center`)/attributes(std_Z)$`scaled:scale`
  # Z_combined = rbind(true_Z, Z)
  # std_Z_combined = scale(Z_combined)
  # std_true_Z = std_Z_combined[1,]
  
  # fitting the true data into our neural netowork...
  pred      <- predict_function(W,std_true_Z) # get the posterior for each parameter
  
  
  # simulation_results[s,] <- obtain_quantiles_and_mean(pred, q_eval) # medians and credible intervals of each posterior
  # unif_rv[s] <- obtain_rank_prob(pred, sim_output$I_summary[choice_of_I])
  if (RESPONSE == 1) this_gammaj = sim_output$I_summary else this_gammaj = sim_output$I_trajectory[round(quantile(1:num_events,RESPONSE))]
  simulation_results_row = c(obtain_quantiles_and_mean(pred, q_eval), 
                             obtain_rank_prob(pred, this_gammaj),
                             this_gammaj)
  simulation_results_row
}
stopCluster(cl)


par(mfrow = c(2,2))

# compare to priors:
true_I_summary = simulation_results[,6]
if (RESPONSE == 1) training_set_hist = I_summaries[[1]] else training_set_hist = I_summaries[[2]][,round(quantile(1:num_events,RESPONSE))]
hist(training_set_hist); hist(true_I_summary, col = rgb(1, 0, 0, 0.5), add = TRUE)
# hist(I_summaries[[2]]); hist(true_I_summaries[[2]], col = rgb(1, 0, 0, 0.5), add = TRUE)
# hist(I_summaries[[3]]); hist(true_I_summaries[[3]], col = rgb(1, 0, 0, 0.5), add = TRUE)
# hist(I_summaries[[4]]); hist(true_I_summaries[[4]], col = rgb(1, 0, 0, 0.5), add = TRUE)

# PIT plot
unif_rv = simulation_results[,5]
qunif.eval = seq(0, 1, by=0.05) # ideal CDF of uniform(0,1) distribution
q.pred = quantile(unif_rv, probs = qunif.eval) # does our distribution have the same quantiles? 
q.unif = qunif.eval

plot(q.unif, q.pred,  type="l", xlab="Theoretical quantile", ylab="Observed quantile",
     cex.lab=1.4, cex.axis=1.4, lwd=2, lty=2, col=choice_of_I+1)
abline(0,1, lty=1, lwd=2)
legend("bottomright", c(expression(I[1]),  expression(I[2]),
                        expression(I[3]), expression(I[4])),
       col = 2:5, lwd=2, lty=2:5)


# compare point estimates with true
plot(simulation_results[,2], true_I_summary, main = "median") # median vs. true 
abline(a=0,b=1)
plot(simulation_results[,4], true_I_summary, main = "mean") # mean vs. true
abline(a=0,b=1)
par(mfrow = c(1,1))

lm_obj_median = lm(true_I_summary~simulation_results[,2])
lm_obj_mean = lm(true_I_summary~simulation_results[,4])

print(sd(true_I_summary))
print(sd(lm_obj_median$residuals))
print(sd(lm_obj_mean$residuals))

R2_median = var(lm_obj_median$fitted.values)/(var(lm_obj_median$fitted.values)+var(lm_obj_median$residuals))
cor_median = sqrt(R2_median)
R2_mean = var(lm_obj_mean$fitted.values)/(var(lm_obj_mean$fitted.values)+var(lm_obj_mean$residuals))
cor_mean = sqrt(R2_mean)

cat("The median estimator explains ", round(R2_median*100,1),"% of the variance.\n")
cat("The mean estimator explains ", round(R2_mean*100,1),"% of the variance.\n")

cat("The correlation between our posterior median and the truth is ", round(cor_median,3),".\n")
cat("The correlation between our posterior mean and the truth is ", round(cor_mean,3),".\n")

coverage = mean((simulation_results[,1] < true_I_summary)&(true_I_summary < simulation_results[,3]))
cat("The overall coverage is ", round(coverage*100,1),".\n")

plot(x = true_I_summary,
     y = simulation_results[,2],
     xlab = paste0("True Infected at time ",round(quantile(1:num_events,.3))),
     ylab = "Posterior Median",
     cex = 1.5,
     cex.lab = 1.6,
     cex.axis = 1.5)
abline(a = 0, b = 1, lwd = 2, col = 2)
points(x = true_I_summary,
       y = simulation_results[,2],cex = 1.5)
# for pic 1... c(100,20)
# 90% credibel interval coverage is 94.6
# R2 is 84.4%
# correlation is .92

# for pic 2... c(120)
# 90% credibel interval coverage is 96.6
# R2 is 89.5%
# correlation is .946

# for pic 3... c(20)
# 90% credibel interval coverage is 89.4
# R2 is 79.4%
# correlation is .891

# for pic 4... c(70)... but need to run twice? 
# 90% credibel interval coverage is 93.4
# R2 is 82.9%
# correlation is .91
