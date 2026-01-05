
# demo inputs ####


# f1 = function(Y){
#   return(Y[1])
# }
# 
# f2 = function(Y){
#   return(Y[2])
# }
# 
# N = 30000
# Y = cbind(rnorm(N, mean = 5),
#           rnorm(N, mean = 5))
# 
# f1_true = apply(Y, 1, f1)
# f2_true = apply(Y, 1, f2)
# 
# theta = rnorm(N, mean = f1_true, sd = f2_true)


# functions ####


act  <- function(x){pmax(x,0)}
actp <- function(x){x>0}

init_hetnorm <- function(p,L1,L2,init_mn=0,init_sd=1,sigma=0.01){
 w       <- list()
 w[[1]]  <- sigma*rnorm(L1)
 w[[2]]  <- sigma*matrix(rnorm(p*L1),p,L1)
 w[[3]]  <- sigma*rnorm(L2)
 w[[4]]  <- sigma*matrix(rnorm(L1*L2),L1,L2)
 w[[5]]  <- init_mn + sigma*rnorm(1)
 w[[6]]  <- sigma*rnorm(L2)
 w[[7]]  <- sigma*rnorm(L1)
 w[[8]]  <- sigma*matrix(rnorm(p*L1),p,L1)
 w[[9]]  <- sigma*rnorm(L2)
 w[[10]] <- sigma*matrix(rnorm(L1*L2),L1,L2)
 w[[11]] <- log(init_sd) + sigma*rnorm(1)
 w[[12]] <- sigma*rnorm(L2)
return(w)}

# test_init_hetnorm = init_hetnorm(p = ncol(Y),
#                                  L1 = 10,
#                                  L2 = 5, 
#                                  init_mn = mean(theta),
#                                  init_sd = sd(theta),
#                                  sigma = .0001)
# length(test_init_hetnorm)
# str(test_init_hetnorm)
# mean(theta)
# mean(test_init_hetnorm[[5]])
# log(sd(theta))
# mean(test_init_hetnorm[[11]])

init_hetnorm2 <- function(p,hidden_layers,init_mn=0,init_sd=1,sigma=0.01){
  node_sequence = c(p,hidden_layers,2) # 2 outputs: mu and sigma
  # initialize with rnorm weights, going layer by layer; sigma is the sd of the rnorm weights
  num_weight_layers = length(node_sequence) - 1
  W = list()
  bias_index = 1 # index for the weight object W[[.]]
  for (layer_index in 1:num_weight_layers){ # allocates the correct number of nodes for each layer
    W[[bias_index]] = sigma*matrix(rnorm(node_sequence[layer_index + 1]), ncol = 1)
    weight_index = bias_index + 1
    W[[weight_index]] = sigma*matrix(rnorm(node_sequence[layer_index]*node_sequence[layer_index + 1]), nrow = node_sequence[layer_index], ncol = node_sequence[layer_index+1])
    bias_index = bias_index + 2
  }
  
  # print(bias_index-2)
  
  W[[bias_index - 2]][1,] = W[[bias_index - 2]][1,] + init_mn
  W[[bias_index - 2]][2,] = W[[bias_index - 2]][2,] + log(init_sd)
  
return(W)}

# test_init_hetnorm2 = init_hetnorm2(p = ncol(Y),
#                                    hidden_layers = c(10,5), 
#                                    init_mn = mean(theta),
#                                    init_sd = sd(theta),
#                                    sigma = .0001)
# length(test_init_hetnorm2)
# str(test_init_hetnorm2)
# mean(theta)
# mean(test_init_hetnorm2[[5]][1,])
# log(sd(theta))
# mean(test_init_hetnorm2[[5]][2,])


predict_hetnorm <- function(w,x){
   m <- act(sweep(x%*%w[[2]],2,w[[1]],"+"))
   m <- act(sweep(m%*%w[[4]],2,w[[3]],"+"))
   m <- m%*%w[[6]] + w[[5]]
   s <- act(sweep(x%*%w[[8]],2,w[[7]],"+"))
   s <- act(sweep(s%*%w[[10]],2,w[[9]],"+"))
   s <- s%*%w[[12]] + w[[11]] 
return(list(mu=m,sigma=exp(s)))}

# test_predict_hetnorm = predict_hetnorm(w = test_init_hetnorm,
#                                        x = Y)
# dim(Y)
# str(test_predict_hetnorm)
# test_predict_hetnorm = predict_hetnorm(w = test_init_hetnorm,
#                                        x = Y[1,])
# str(test_predict_hetnorm)

predict_hetnorm2 <- function(w,x){
  
  # forward pass the inputs X through the hidden layers (linear combo, bias, activation function)
  for (i in seq(1, length(w)-2, by = 2)){
    x = sweep(x%*%w[[i+1]], 2, w[[i]], "+")
    x = act(x)
  }
  
  # one last forward pass, but don't apply the activation function
  next_index = length(w) - 1
  x = sweep(x%*%w[[next_index+1]] , 2, w[[next_index]], "+")
  
return(list(mu=x[,1],sigma=exp(x[,2])))}

# test_predict_hetnorm2 = predict_hetnorm2(w = test_init_hetnorm2,
#                                          x = Y)
# dim(Y)
# str(test_predict_hetnorm2)
# test_predict_hetnorm2 = predict_hetnorm2(w = test_init_hetnorm2,
#                                          x = Y[1,])
# str(test_predict_hetnorm2)


loss_hetnorm <- function(w,x,y){
   l <- predict_hetnorm(w,x)
   l <- -sum(dnorm(y,l$mu,l$sigma,log=TRUE))
return(l)}

# test_loss_hetnorm = loss_hetnorm(w = test_init_hetnorm,
#                                  x = Y,
#                                  y = theta)
# test_loss_hetnorm

loss_hetnorm2 <- function(w,x,y){
  l <- predict_hetnorm2(w,x)
  l <- -sum(dnorm(y,l$mu,l$sigma,log=TRUE))
  return(l)}

# test_loss_hetnorm2 = loss_hetnorm2(w = test_init_hetnorm2,
#                                    x = Y,
#                                    y = theta)
# test_loss_hetnorm2


grad_hetnorm2 <- function(w,x,y,transfer_learning = FALSE,priorweights = rep(1,length(y))){
  
  # forward pass ####
  
  previous_nodes = list() # previous_node[[layer2]] is node[[layer1]] 
  previous_nodes[[1]] = x # 
  unactivated_nodes = list() # each layer's unactivated_nodes are used in the Backpropagation derivative
  for (l in seq(1, length(w)-2, by = 2)){ # forward pass through the NN's hidden layers...
    current_layer = ceiling(l/2)
    # print(dim(x))
    # print(dim(w[[l+1]]))
    unactivated_nodes[[current_layer]] = sweep(x%*%w[[l+1]], 2, w[[l]], "+") # save
    x = sweep(x%*%w[[l+1]], 2, w[[l]], "+") # update X and continue to activation
    previous_nodes[[current_layer + 1]] = act(x) # save
    x = act(x) # update X and continue to the next layer
  } # ...after this for-loop, do the last output layer
  l = l + 2
  current_layer = ceiling(l/2)
  unactivated_nodes[[current_layer]] = 1 # output layer's chain rule is "1" in Backpropagation
  output = sweep(x%*%w[[l+1]] , 2, w[[l]], "+") # save this for Stage 2
  
  
  # derivatives of the loss ####
  
  # print(dim(output))
  
  m3 <- output[,1]
  s3 <- output[,2]
  
  dldm3  <- -exp(-2*s3)*(y-m3)
  dlds3  <- 1+dldm3*(y-m3)
  
  dldoutput <- cbind(dldm3, dlds3)
  dldoutput <- sweep(dldoutput, 1, priorweights, "*") # multiply by the VaNBayes weights
  
  # print(dim(dldoutput))
  
  # backpropagation ####
  
  # our calculated gradient will be the same dimension as the inputted weights
  g  <- w
  
  # this is required to fit the for-loop
  w[[length(w) + 2]] = diag(2)
  
  # the backpropagation algorithm
  # more details in the Overleaf document
  # Idea: Uses 
  # Requirements: 
  # 1) the previous layer's nodes (the coefficients of W, so when the derivative of WX is taken wrt W, we need X)
  # 2) the current node's unactivated layers (when we take the derivative of act(unactivated), we calculate act_derivative(unactivated))
  # 3) the previously considered weights (since they are coefficients for the part we are taking the derivative of)
  # 4) the dlossdoutput, which is the beginning of the chain rule product
  # we continuously update the "left-side" Nxp factor as we go through layers (chain rule)
  
  last_layer_to_train = ifelse(transfer_learning, length(w)-2, 1)
  
  left_side <- dldoutput # start as dldoutput
  for (l in seq(length(w)-2, last_layer_to_train, by = -2)){ # notice we go in opposite direction as when going through the prediction
    current_layer = ceiling(l/2)
    # cat("the current trained layer is", current_layer, "\n")
    left_side <- (left_side%*%t(w[[l+2]]))*actp(unactivated_nodes[[current_layer]]) # update left-side factor from Nx(p1) to Nx(p2)
    g[[l]] <- t(t(left_side) %*% previous_nodes[[current_layer]]) # the previous node as a coefficient
    g[[l-1]] <- t(left_side) %*% rep(1, ncol(t(left_side))) # no coefficient for the bias's derivative, so just "1" 
  }
  
  return(g) # we are finally done!
  
}

# test_grad_hetnorm2 = grad_hetnorm2(w = test_init_hetnorm2,
#                                    x = Y,
#                                    y = theta)
# test_grad_hetnorm2




if(FALSE){ # Check gradient

 p  <- 3
 L1 <- 5
 L2 <- 2
 n  <- 7
 x  <- matrix(rnorm(n*p),n,p)
 y  <- rnorm(n)

 w   <- init_hetnorm2(p,c(L1,L2))
 l0  <- loss_hetnorm2(w,x,y)

 numgrad <- function(w,x,y,k,l0,eps){
   g <- w[[k]]
   if(is.vector(g)){
      for(i in 1:length(g)){
        We         <- w
        We[[k]][i] <- We[[k]][i] + eps
        g[i]       <- loss_hetnorm2(We,x,y)
      }
   }
   if(is.matrix(g)){
      for(i in 1:nrow(g)){for(j in 1:ncol(g)){
        We <- w
        We[[k]][i,j] <- We[[k]][i,j] + eps
        g[i,j] <- loss_hetnorm2(We,x,y)
      }}
   }
 return((g-l0)/eps)}

 
 eps <- 10^(-8)
 g1   <- w
 for(k in 1:length(g1)){g1[[k]]<-numgrad(w,x,y,k,l0,eps)}
 g2 <- grad_hetnorm2(w,x,y)
 print(g1[[7]])
 print(g2[[7]])
 for(k in 1:length(g1)){
   print(k)
   print(g1[[k]])
   print(g2[[k]])
 }
}

