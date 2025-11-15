# sim_study_marginal_PLOTS


# SIMULATION STUDY ####

## Simulation Settings ####
# settings
true_lambda= .4
true_I0 = 20
true_psi = 5
true_parameters = c(true_lambda, true_I0, true_psi)
parameter_symbols = c(expression(lambda), expression(I[0]), expression(psi))

# number of simulation study datasets
nsims = 100

# which quantiles to compute
q_eval = c(0.05, 0.5, 0.95) 

## Simulation Study ####
simulation_results = array(dim = c(nsims, 3, length(q_eval), length(N_vec)))
sim_study_datasets = matrix(nrow = nsims, ncol = length(configure_data(run_siminf_fast(true_lambda, true_I0, true_psi))))
for(k in 1:length(N_vec)){
  N = N_vec[k]
  print(N)
  # TRAINING ####
  
  ## Generate nruns of data used for VaNBayes ####
  size_N_dataset = generate_Z(N)
  Z = size_N_dataset$data
  theta = size_N_dataset$parameters
  
  ### configure ####
  gamma = theta; for (i in 1:nrow(theta)) gamma[i,] = configure_parameters(theta[i,])
  
  ## Fit VaNBayes to the generated data ####
  
  # fit the neural networks (one for each parameter)
  # response for the jth neural network
  gammaj     <- gamma[,choice_of_theta] 
  
  # estimate the marginals using the normal variational posterior (assume heterogeneous normal)
  if(which_arch == 1){
    initj  <- init_function(p=ncol(Z),
                            hidden_layers = architecture,
                            init1,
                            init2)
    
  } else {
    initj      <- init_function(p=ncol(Z),
                               L1=architecture[1],
                               L2=architecture[2],
                               init1,
                               init2) # initialize weights with normal distribution
  }
  
  # adaptive SGD using our standardized PC Scores and the log-normal as the loss function
  modelj     <- adam(w=initj, 
                     x=Z, 
                     y=gammaj, 
                     loss=loss_function, 
                     grad=grad_function,
                     batchsize = BATCHSIZE[k], 
                     epochs = EPOCHS, 
                     propval = PROPVAL,
                     early_stop = EARLY_STOP,
                     lr = LR,
                     verbose = 2)
  
  # save the weights of each neural network
  W = as.list(modelj$w)
  
  for(sim in 1:nsims){
    print(sim)
    set.seed(sim)
    Y0 <- run_siminf_fast(true_lambda, true_I0, true_psi)
    Z0 <- configure_data(Y0)
    sim_study_datasets[sim,] <- Z0
    
    # predict
    pred      <- predict_function(W,Z0) # get the posterior for each parameter
    if (var_post == "gamma"){
      simulation_results[sim,choice_of_theta,,k] <- qgamma(q_eval,shape = pred$a,scale = pred$b) # medians and credible intervals of each posterior
    } else {
      simulation_results[sim,choice_of_theta,,k] <- qnorm(q_eval,mean = pred$mu,sd = pred$sigma) # medians and credible intervals of each posterior
    }
  }
  
}
write.csv(sim_study_datasets, file = "./sim_study_datasets.csv")



## Transforming and Saving VaNBayes results ####
true_gamma = configure_parameters(true_parameters)

# coverage and MAD
parameter_coverage = matrix(nrow = length(true_parameters), ncol = length(N_vec))
parameter_MAD = matrix(nrow = length(true_parameters), ncol = length(N_vec))
parameter_SDAD = matrix(nrow = length(true_parameters), ncol = length(N_vec))
for (k in 1:length(N_vec)){
  # coverage
  parameter_coverage[,k] = c(sum((simulation_results[,1,1,k] < true_gamma[1]) & (true_gamma[1] < simulation_results[,1,3,k]))/nsims,
                             sum((simulation_results[,2,1,k] < true_gamma[2]) & (true_gamma[2] < simulation_results[,2,3,k]))/nsims,
                             sum((simulation_results[,3,1,k] < true_gamma[3]) & (true_gamma[3] < simulation_results[,3,3,k]))/nsims)
  # MAD
  parameter_MAD[,k] = c(mean(abs(simulation_results[,1,2,k] - true_gamma[1])), 
                        mean(abs(simulation_results[,2,2,k] - true_gamma[2])), 
                        mean(abs(simulation_results[,3,2,k] - true_gamma[3])))
  
  # SD
  parameter_SDAD[,k] = c(sd(abs(simulation_results[,1,2,k] - true_gamma[1])), 
                         sd(abs(simulation_results[,2,2,k] - true_gamma[2])), 
                         sd(abs(simulation_results[,3,2,k] - true_gamma[3])))
  
  
}

print("VaNBayes Parameter Coverage:")
print(parameter_coverage[choice_of_theta,])
print("VaNBayes Parameter MAD")
print(parameter_MAD[choice_of_theta,])
print("VaNBayes Parameter SDAD")
print(parameter_SDAD[choice_of_theta,])

VaNBayes_medians = simulation_results[,,2,]


spacing = 0.1
plot(x = 1:length(N_vec),
     y = rep(NA,length(N_vec)),
     axes = F,
     cex.lab=1.5,
     main = parameter_symbols[choice_of_theta],
     xlab = "Number of Simulated Samples",
     ylab = "Posterior Medians",
     xlim = range(c(1:length(N_vec)-0.25, 1:length(N_vec)+0.25)),
     # xlim = range(c(min(N_vec)-0.25, max(N_vec)+0.25)),
     ylim = range(c(VaNBayes_medians[,choice_of_theta,],true_gamma[choice_of_theta])))
abline(h = 0, col = "black", lwd = 2)
for (k in 1:length(N_vec)){
  boxplot(VaNBayes_medians[,choice_of_theta,k],at=k-spacing,add=TRUE,axes=F,col="lightblue",boxwex=.4,outline=FALSE)
}
abline(h = true_gamma[choice_of_theta], col = "red")
axis(2,cex.axis=1.25,cex.lab=1.25)
axis(1,at=1:length(N_vec),lab=N_vec/2,cex.axis=1.5,cex.lab=1.5)
# axis(1,at=N_vec,lab=N_vec,cex.axis=1.5,cex.lab=1.5)
legend("topright", inset=c(0,-.1),xpd=TRUE, # need to adjust the first inset number
       legend = c("VaNBayes"),
       fill=c("lightblue"),
       bty="n",cex=1.5)

