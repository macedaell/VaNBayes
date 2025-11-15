
script_to_test = "NegativeBinomial1.R"

source(script_to_test)
source("../adam.R")

# results: I don't think the third parametrization works well...

if (script_to_test == "NegativeBinomial1.R"){
  
  N <- 10000
  theta <- rnbinom(10000, prob = 0.5, size = 10)
  theta_val <- rnbinom(10000, prob = 0.5, size = 10)
  
  
  n <- 5
  
  Y <- matrix(nrow = N, ncol = n)
  Y_val <- matrix(nrow = N, ncol = n)
  for (i in 1:nrow(Y)){
    Y[i,] = rnorm(n, mean = theta[i], sd = 1)
    Y_val[i,] = rnorm(n, mean = theta_val[i], sd = 1)
  }
  
  
  MOM_r = mean(theta)/(var(theta) - mean(theta))
  MOM_p = mean(theta)/var(theta)
  
  init <- init_NegativeBinomial1(p = ncol(Y),
                                hidden_layers = c(5),
                                init_r = MOM_r,
                                init_p = MOM_p)
  
  loss_function = loss_NegativeBinomial1
  grad_function = grad_NegativeBinomial1
  predict_function = predict_NegativeBinomial1
  
  model <- adam(w = init,
                x = Y,
                y = theta, 
                loss = loss_function,
                grad = grad_function,
                batchsize = 100,
                epochs = 200,
                propval = 0.2, 
                early_stop = 5,
                lr = 0.001,
                verbose = 2)
  W <- as.list(model$w)
  
  
  pred <- predict_function(W, Y_val)
  
  final_prediction = qnbinom(0.5,size = pred$r, prob = pred$p)
  
  
} else if (script_to_test == "NegativeBinomial2.R"){
  
  N <- 10000
  theta <- rnbinom(10000, mu = 10, size = 10)
  theta_val <- rnbinom(10000, mu = 10, size = 10)
  
  
  n <- 5
  
  Y <- matrix(nrow = N, ncol = n)
  Y_val <- matrix(nrow = N, ncol = n)
  for (i in 1:nrow(Y)){
    Y[i,] = rnorm(n, mean = theta[i], sd = 1)
    Y_val[i,] = rnorm(n, mean = theta_val[i], sd = 1)
  }
  
  
  init <- init_NegativeBinomial2(p = ncol(Y),
                                hidden_layers = c(5),
                                init_mn = mean(theta),
                                init_sd = sd(theta))
  
  loss_function = loss_NegativeBinomial2
  grad_function = grad_NegativeBinomial2
  predict_function = predict_NegativeBinomial2
  
  model <- adam(w = init,
                x = Y,
                y = theta, 
                loss = loss_function,
                grad = grad_function,
                batchsize = 100,
                epochs = 200,
                propval = 0.2, 
                early_stop = 5,
                lr = 0.001,
                verbose = 2)
  W <- as.list(model$w)
  
  
  pred <- predict_function(W, Y_val)
  
  final_prediction = pred$mu
  
} else {
  
  N <- 10000
  theta <- rnbinom(10000, mu = 10, size = 10)
  theta_val <- rnbinom(10000, mu = 10, size = 10)
  
  
  n <- 5
  
  Y <- matrix(nrow = N, ncol = n)
  Y_val <- matrix(nrow = N, ncol = n)
  for (i in 1:nrow(Y)){
    Y[i,] = rnorm(n, mean = theta[i], sd = 1)
    Y_val[i,] = rnorm(n, mean = theta_val[i], sd = 1)
  }
  
  init <- init_NegativeBinomial3(p = ncol(Y),
                                 hidden_layers = c(10,5),
                                 init_mn = mean(theta))
  
  loss_function = loss_NegativeBinomial3
  grad_function = grad_NegativeBinomial3
  predict_function = predict_NegativeBinomial3
  
  model <- adam(w = init,
                x = Y,
                y = theta, 
                loss = loss_function,
                grad = grad_function,
                batchsize = 100,
                epochs = 200,
                propval = 0.2, 
                early_stop = 5,
                lr = 0.001,
                verbose = 2)
  W <- as.list(model$w)
  
  
  pred <- predict_function(W, Y_val)
  
  final_prediction = pred$mu
  
}


plot(theta_val, final_prediction)
abline(a=0,b=1)







