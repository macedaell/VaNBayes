
# SETUP ####

# training arguments
N_vec = c(500, 1000, 1500, 2000, 2500, 3000, 3500, 4000, 4500, 5000, 5500, 6000) # Brian liked something like: c(100, 500, 1000, 2000, 3000, 4000, 5000, 6000)
PROPVAL = 0.5 # 0.2

# VaNBayes arguments
BATCHSIZE = rep(100,length(N_vec)) # next, make this bigger
EPOCHS = 200
EARLY_STOP = 30

set.seed(1)

source("normal.R")
source("gamma2.R")
source("adam.R")
source("hetnorm2.R")
source("aux_functions.R")
library(SimInf)
library(Matrix)
library(truncnorm)

library(doParallel)
options(warn = 2)


original_parameters = function(configured_parameters){
  parameters = configured_parameters
  parameters[1] = boxcox_inv(configured_parameters[1],0.1)
  parameters[2] = boxcox_inv(configured_parameters[2],0.25)
  # parameters[3] = boxcox_inv(configured_parameters[3],0.15)
  parameters[3] = exp(configured_parameters[3]) - 1 
  return(parameters)
}

configure_parameters = function(parameters){
  configured_parameters = parameters
  # # OPTION 1: change to standard normal distribution
  # configured_parameters[1] = (log(parameters[1]) - log(.4))/.5
  # configured_parameters[2] = qnorm(pgamma(parameters[2], 2, scale = 20))
  # configured_parameters[3] = qnorm(pexp(parameters[3], 1/5))
  
  # # OPTION 2: small transformations
  # configured_parameters[1] = log(parameters[1])
  configured_parameters[1] = boxcox(parameters[1],0.1)
  # configured_parameters[2] = log(parameters[2]) # roughly normal
  configured_parameters[2] = boxcox(parameters[2],0.25) # roughly normal
  # configured_parameters[2] = (boxcox(parameters[2],-0.7)) # roughly normal
  # configured_parameters[2] = qt(pgamma(parameters[2], 2, scale = 20),5) # OK, not great but seems fixable
  # configured_parameters[2] = parameters[2]
  # configured_parameters[3] = boxcox(parameters[3],0.15)
  configured_parameters[3] = log(parameters[3] + 1) # this is a little better(?) if we want to show VaNBayes doing worse
  return(configured_parameters)
}

boxcox = function(x,lambda) (x^lambda - 1)/lambda
boxcox_inv = function(y,lambda) (lambda*y + 1)^(1/lambda)

# FITTING THE NEURAL NETWORKS IN GENERAL ####

choice_of_theta = 3 # theta3 is okay. Could be better...
# do this same procedure with all thetas until the results look good enough to share with Brian
# maybe sometime today ask Brian what kinds of plots to use? 

N = N_vec[1]

var_posts =             c("normal",        "normal",  "normal")
architectures = list(c(10, 5), c(20, 10), c(25,15)) # 25,15 is pretty good
LRs = c(0.001, 0.001, 0.001) # 0.001 is fine... 0.0004 is high quality

LR = LRs[choice_of_theta]
var_post = var_posts[choice_of_theta]
architecture = architectures[[choice_of_theta]]


size_N_dataset = generate_Z(N)

theta = size_N_dataset$parameters

Z = size_N_dataset$data

gamma = theta; for(i in 1:nrow(theta)) gamma[i,] = configure_parameters(theta[i,])

plot(density(gamma[,choice_of_theta]))

if (var_post == "gamma"){
  init_function = init_gamma2
  init1 = 1.185  # might need to customize these
  init2 = 1/0.95 # 
  loss_function = loss_gamma2
  grad_function = grad_gamma2
  predict_function = predict_gamma2
  which_arch = 1
} else if (var_post == "normal"){
  init_function = init_hetnorm2
  init1 = mean(gamma[,choice_of_theta])
  init2 = sd(gamma[,choice_of_theta])
  loss_function = loss_hetnorm2
  grad_function = grad_hetnorm2
  predict_function = predict_hetnorm2
  which_arch = 1
} else { # hetnorm
  init_function = init_hetnorm
  init1 = mean(gamma[,choice_of_theta])
  init2 = sd(gamma[,choice_of_theta])
  loss_function = loss_hetnorm
  grad_function = grad_hetnorm
  predict_function = predict_hetnorm
  which_arch = 2
}

# Estimate the weights
gammaj <- gamma[,choice_of_theta]
if(which_arch == 1){
  initj  <- init_function(p=ncol(Z),
                          hidden_layers = architecture,
                          init1,
                          init2)
  
} else { # only when var_post == "hetnorm"
  initj  <- init_function(p=ncol(Z),
                          L1 = architecture[1], L2 = architecture[2],
                          init1,
                          init2)
  
}
modelj <- adam(w=initj, 
               x=Z, 
               y=gammaj,
               loss=loss_function, 
               grad=grad_function,
               epochs = EPOCHS, 
               batchsize = BATCHSIZE[1],
               early_stop = EARLY_STOP,
               lr = LR,
               propval = PROPVAL,
               verbose = 2)
W = as.list(modelj$w)


# PIT PLOTS ####



val_dataset = generate_Z(10000)
theta_val = val_dataset$parameters
Z_val = val_dataset$data

gamma_val = theta_val; for(i in 1:nrow(theta_val)) gamma_val[i,] = configure_parameters(theta_val[i,])

unif_rv = matrix(nrow = nrow(Z_val), ncol = 3)

for (k in 1:nrow(Z_val)){
  
  pred      <- predict_function(W,Z_val[k,])
  if (var_post == "gamma"){
    unif_rv[k,choice_of_theta] <- pgamma(gamma_val[k,choice_of_theta],shape = pred$a,scale = pred$b)
  } else {
    unif_rv[k,choice_of_theta] <- pnorm(gamma_val[k,choice_of_theta],mean = pred$mu,sd = pred$sigma)
  }
}


par(mfrow = c(1,2))
hist(unif_rv[,choice_of_theta]) # should look uniform...

qunif.eval = seq(0, 1, by=0.05) # ideal CDF of uniform(0,1) distribution
q.pred = quantile(unif_rv[,choice_of_theta], probs = qunif.eval) # does our distribution have the same quantiles? 
q.unif = qunif.eval

plot(q.unif, q.pred,  type="l", xlab="Theoretical quantile", ylab="Observed quantile",
     cex.lab=1.4, cex.axis=1.4, lwd=2, lty=2, col=choice_of_theta+1)
abline(0,1, lty=1, lwd=2)
legend("bottomright", c(expression(theta[1]),  expression(theta[2]),
                        expression(theta[3])),
       col = 2:4, lwd=2, lty=2:4)
par(mfrow = c(1,1))

# source("sim_study_marginal.R")
# source("sim_study_all.R")
