# used to simulate observed data
negative_binom_fast <- function(y, dispersion){
  
  r = dispersion
  var = y + (1/r) * y^2
  p = (var-y)/var
  return(rnbinom(length(p), size = r, prob = 1-p))
  
}

# add the np.clip function, to be used in the SIR simulation
numpy_clip = function(x, a, b){
  return(ifelse(x <= a, a, ifelse(x>=b,b,x)))
}

# simulate from the SIR model
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
  # simdata = log(simdata+1)%*%loading_matrix
  return(simdata)
  
}

# used to make a new dataset every time a new sample size is considered
generate_data = function(nruns){
  demo_y <- run_siminf_fast(0.4, 20, 5)
  generated_data = matrix(0, nrow = nruns, ncol = length(demo_y)+3) # 14 infected days, 3 parameters
  for (i in 1:nruns) {
    # if (i %% 10 == 0) cat(paste("\rOn iter", i," out of ", nruns))
    
    # draw from priors
    my_lambd <- rlnorm(1, log(.4), .5)
    my_I0 <- rgamma(1,2,scale = 20)
    my_psi <- rexp(1, 1/5)
    
    # generate an SIR dataset
    cand_y <- run_siminf_fast(my_lambd, my_I0, my_psi)
    
    # save the results
    generated_data[i,] = c(my_lambd,my_I0,my_psi,c(cand_y))
  }
  cat("\n")
  
  parameters = generated_data[,1:3]
  data = generated_data[,-c(1:3)]
  
  results = list(parameters = parameters,
                 data = data)
  
}



# PCA LOADING STUFF ####

# CREATE THE LOADING MATRIX FOR Y ####
size_N_dataset = generate_data(10000)
Y = size_N_dataset$data
scaledlogY = scale(log(Y+1))
logY_mn = attributes(scaledlogY)$`scaled:center`
logY_sd = attributes(scaledlogY)$`scaled:scale`
Loading_scaledlogY = eigen(cov(scaledlogY))$vectors[,1:7] # 90% var. expl.
scaledZ = scale(scaledlogY%*%Loading_scaledlogY)
plot(density(scaledZ)) # looks good! Nice and symmetric. No crazy values
Z_mn = attributes(scaledZ)$`scaled:center`
Z_sd = attributes(scaledZ)$`scaled:scale`

write.csv(Loading_scaledlogY, file = "./loading_matrix.csv")
write.csv(logY_mn, file = "./logY_mn.csv")
write.csv(logY_sd, file = "./logY_sd.csv")
write.csv(Z_mn, file = "./Z_mn.csv")
write.csv(Z_sd, file = "./Z_sd.csv")


loading_matrix = as.matrix(read.csv("./loading_matrix.csv", )[,-1])
colnames(loading_matrix) = NULL


configure_data = function(data){
  scaledlogY = (log(data+1) - logY_mn)/logY_sd
  Z = scaledlogY%*%loading_matrix
  scaledZ = (Z - Z_mn)/Z_sd
  return(scaledZ)
}



generate_Z = function(nruns){
  demo_y <- run_siminf_fast(0.4, 20, 5)
  demo_y <- configure_data(demo_y)
  generated_data = matrix(0, nrow = nruns, ncol = length(demo_y)+3) # 14 infected days, 3 parameters
  for (i in 1:nruns) {
    # if (i %% 10 == 0) cat(paste("\rOn iter", i," out of ", nruns))
    
    # draw from priors
    my_lambd <- rlnorm(1, log(.4), .5)
    my_I0 <- rgamma(1,2,scale = 20)
    my_psi <- rexp(1, 1/5)
    
    # generate an SIR dataset
    cand_y <- run_siminf_fast(my_lambd, my_I0, my_psi)
    cand_z <- configure_data(cand_y)
    
    # save the results
    generated_data[i,] = c(my_lambd,my_I0,my_psi,c(cand_z))
  }
  cat("\n")
  
  parameters = generated_data[,1:3]
  data = generated_data[,-c(1:3)]
  
  results = list(parameters = parameters,
                 data = data)
  
}

Z_test = generate_Z(10000)$data
plot(density(Z_test)) # finally looks good!
