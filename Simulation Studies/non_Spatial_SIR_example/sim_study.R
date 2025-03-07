
# Zika model fit via ABC
# use "zika Brazil.rdata" from main data analysis

rm(list = ls())
source("adam.R")
source("hetnorm2.R")
library(SimInf)
library(Matrix)
library(truncnorm)
options(warn = 2)

#########
# Setup #
#########

nruns <- 100000 # should be around 100,000


negative_binom_fast <- function(y, dispersion){
  
  r = dispersion
  var = y + (1/r) * y^2
  p = (var-y)/var
  return(rnbinom(length(p), size = r, prob = 1-p))
  # return(rnbinom(length(y), size = dispersion, mu = y)) # using mu, either seem to work
}


numpy_clip = function(x, a, b){
  return(ifelse(x <= a, a, ifelse(x>=b,b,x)))
}


run_siminf_fast <- function(my_lambda, I0, dispersion) {
  
  # global variable
  N = 83e6
  recover = .1
  
  # initial conditions
  S = c(N - I0)
  I = c(I0)
  R = c(0)
  
  # reported new cases
  C = c(I0)
  
  for (i in 1:13){ # needs to be 13!
    
    # calculate new cases
    I_new = my_lambda*(I[length(I)]*S[length(S)]/N)
    
    # SIR equations
    S_t = S[length(S)] - I_new
    I_t = numpy_clip(I[length(I)] + I_new - recover*I[length(I)], 0, N)
    R_t = numpy_clip(R[length(R)] + recover*I[length(I)], 0, N)
    
    S = c(S, S_t)
    I = c(I, I_t)
    R = c(R, R_t)
    C = c(C, I_new)
  }
  # print(C)
  # print(dispersion)
  simdata = negative_binom_fast(numpy_clip(C, 0, N) + 1e-5,dispersion)
  
  return(simdata)
  
}



##################################################################################
# Generating The Data from the Training Distribution to train the NN in VaNBayes #
##################################################################################

results = matrix(0, nrow = nruns, ncol = 14+3) # 14 infected days, 3 parameters

for (i in 1:nruns) {
  if (i %% 10 == 0) print(paste("On iter", i))
  
  # draw from priors
  my_lambd <- rlnorm(1, log(.4), .5)
  my_I0 <- rgamma(1,2,scale = 20)
  my_psi <- rexp(1, 1/5)
  
  # generate an SIR dataset
  cand_y <- run_siminf_fast(my_lambd, my_I0, my_psi)

  # save the results
  results[i,] = c(my_lambd,my_I0,my_psi,c(cand_y))
}

parameters = results[,1:3]
data = results[,-c(1:3)]



######################### FITTING THE NEURAL NETWORKS ##########################

# transform the parameters so the heterogeneous normal model makes sense for them


# lambda--made from lognormal
normal_lambda = (log(parameters[,1]) - log(.4))/.5

# I0-- made from gamma
uniform_I0 = pgamma(parameters[,2], 2, scale = 20)
normal_I0 = qnorm(uniform_I0)

# psi-- made from exponential
uniform_psi = pexp(parameters[,3], 1/5)
normal_psi = qnorm(uniform_psi)

normal_parameters = cbind(normal_lambda, normal_I0, normal_psi)

scaled_data = scale(data)

# Neural Network Architecture -- completely connected with two layers
layer1_nodes = 200
layer2_nodes = 50

# fit the neural networks (one for each parameter)
W  <- list()
for(j in 1:ncol(normal_parameters)){
  print(paste("Training parameter",j))
  
  gammaj     <- normal_parameters[,j] # response for the jth neural network
  
  
  initj      <- init_hetnorm(p=ncol(scaled_data),L1=layer1_nodes,L2=layer2_nodes,
                             init_mn=mean(gammaj),init_sd=sd(gammaj)) # initialize weights with normal distribution
  
  # adaptive SGD using our standardized PC Scores and the log-normal as the loss function
  modelj     <- adam(w=initj, x=scaled_data, y=gammaj, 
                     loss=loss_hetnorm, grad=grad_hetnorm, verbose = 2) 
  
  W[[j]] <- as.list(modelj$w) # save the weights of each neural network
}


########################### SIMULATION STUDY ###################################

# simulation setting
true_lambda= .4
true_I0 = 20
true_psi = 5



###### change everything below this

true_parameters = c(true_lambda, true_I0, true_psi)

nsims = 100
q_eval = c(0.05, 0.5, 0.95) # medians for the point estimates, and 0.05, 0.95 for the 90% credible intervals

simulation_results = array(0, c(nsims, 3, length(q_eval)))
sd_results = matrix(0, nrow = nsims, ncol = 3)
sim_study_datasets = matrix(nrow = nsims, ncol = ncol(data))
for(sim in 1:nsims){
  print(sim)
  set.seed(sim)
  
  example_dataset <- run_siminf_fast(true_lambda, true_I0, true_psi)
  sim_study_datasets[sim,] <- example_dataset
  
  # fitting the true data into our neural network...
  for (j in 1:length(true_parameters)){
    
    # standardizing the dataset
    standardized_example_df = (example_dataset - attributes(scaled_data)$`scaled:center`)/attributes(scaled_data)$`scaled:scale`
    
    # the VaNBayes prediction
    pred      <- predict_hetnorm(W[[j]],standardized_example_df) # get the posterior for each parameter
    simulation_results[sim,j,] <- qnorm(q_eval,pred$mu,pred$sigma) # medians and credible intervals of each posterior
    sd_results[sim,j] = pred$sigma[1,1] # standard deviation of each posterior 
  }
}
write.csv(sim_study_datasets, file = "sim_study_datasets.csv")


######## Transforming into original scale

# transform lambdas back to prior-scale
simulation_results[,1,] = exp(simulation_results[,1,]*.5 + log(.4))

# transform I0 back to prior-scale
simulation_results[,2,] = qgamma(pnorm(simulation_results[,2,]), 2, scale = 20)

# transform psi back to prior-scale
simulation_results[,3,] = qexp(pnorm(simulation_results[,3,]), 1/5)

# coverage
parameter_coverage = c(sum((simulation_results[,1,1] < true_lambda) & (true_lambda < simulation_results[,1,3]))/nsims,
                              sum((simulation_results[,2,1] < true_I0) & (true_I0 < simulation_results[,2,3]))/nsims,
                              sum((simulation_results[,3,1] < true_psi) & (true_psi < simulation_results[,3,3]))/nsims)

# MAD
parameter_MAD = c(mean(abs(simulation_results[,1,2] - true_lambda)), 
                         mean(abs(simulation_results[,2,2] - true_I0)), 
                         mean(abs(simulation_results[,3,2] - true_psi)))

# boxplot of posterior medians vs. true values
labs <- c(expression(lambda), expression(I[0]), expression(psi))
c1 <- c(parameter_coverage[1], parameter_coverage[2], parameter_coverage[3])
boxplot(simulation_results[,,2],axes=F,outline=F, ylab="Posterior median",cex.lab=1.25)
axis(2,cex.axis=1.25,cex.lab=1.25)
axis(1,at=1:3,lab=labs,cex.axis=1.25,cex.lab=1.25)
points(true_parameters,col=2,pch=19)
text(1:3,-0.05,c1,pos=3,cex=1.25)

# import bayesflow results for comparison
bayesflow_medians = read.csv("bayesflow_posterior_medians.csv", header = FALSE)

# boxplot of bayesflow posterior medians vs. true values
labs <- c(expression(lambda), expression(I[0]), expression(psi))
boxplot(bayesflow_medians,axes=F,outline=F, ylab="Posterior median",cex.lab=1.25)
axis(2,cex.axis=1.25,cex.lab=1.25)
axis(1,at=1:3,lab=labs,cex.axis=1.25,cex.lab=1.25)
points(true_parameters,col=2,pch=19)


# combined boxplot
plot(NA,xlim=c(0,3+1), ylim=range(c(bayesflow_medians, simulation_results[,,2])),axes=F,ylab="Posterior Median",cex.lab=1.5)
for(j in 1:3){ # also might need to make this color contrast stronger...
  boxplot(simulation_results[,j,2],at=j-0.25,add=TRUE,axes=F,col="white",boxwex=.4,outline=FALSE)
  # boxplot(p1[,j],at=j,add=TRUE,axes=F,col= "lightgray",boxwex=.4,outline=FALSE)
  points(x = j, y = true_parameters[j], col=2,pch=19)
  boxplot(bayesflow_medians[,j],at=j+0.25,add=TRUE,axes=F,col= "darkgray",boxwex=.4,outline=FALSE)
}
axis(2,cex.axis=1.25,cex.lab=1.25)
axis(1,at=1:3,lab=labs,cex.axis=1.25,cex.lab=1.25)
legend("topright",c("VaNBayes", "BayesFlow"),fill=c("white","darkgray"),bty="n",cex=1.5)


mytitles = c("infection_rate", "initial_infected_counts", "dispersion_param")
# faceted boxplots
for (j in 1:3){
  boxplot(cbind(simulation_results[,j,2], bayesflow_medians[,j]),axes=F,outline=F, ylab="Posterior median",cex.lab=1.25, main = mytitles[j])
  labs <- c(expression(lambda), expression(I[0]), expression(psi))
  abline(h = true_parameters[j], col = "red", lwd = 2)
  # boxplot(bayesflow_medians,axes=F,outline=F, ylab="Posterior median",cex.lab=1.25)
  axis(2,cex.axis=1.25,cex.lab=1.25)
  axis(1,at=1:2,lab=c("VaNBayes", "Bayesflow"),cex.axis=1.25,cex.lab=1.25)
}


