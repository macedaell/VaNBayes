
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

logit <- function(x){log(x)-log(1-x)}
expit <- function(x){1/(1+exp(-x))}
act  <- function(x){pmax(x,0)}
actp <- function(x){x>0}


init_NegativeBinomial1 <- function(p,hidden_layers,init_r=0,init_p=1,sigma=0.01){
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
  
  W[[bias_index - 2]][1,] = W[[bias_index - 2]][1,] + init_r
  W[[bias_index - 2]][2,] = W[[bias_index - 2]][2,] + init_p
  
return(W)}


predict_NegativeBinomial1 <- function(w,x){
  
  # forward pass the inputs X through the hidden layers (linear combo, bias, activation function)
  for (i in seq(1, length(w)-2, by = 2)){
    x = sweep(x%*%w[[i+1]], 2, w[[i]], "+")
    x = act(x)
  }
  
  # one last forward pass, but don't apply the activation function
  next_index = length(w) - 1
  x = sweep(x%*%w[[next_index+1]] , 2, w[[next_index]], "+")
  
return(list(r=exp(x[,1]),p=expit(x[,2])))}


loss_NegativeBinomial1 <- function(w,x,y){
  l <- predict_NegativeBinomial1(w,x)
  l <- -sum(dnbinom(y,size=l$r,prob=l$p,log=TRUE))
  return(l)}


grad_NegativeBinomial1 <- function(w,x,y,transfer_learning = FALSE,priorweights = rep(1,length(y))){
  
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
  
  dldm3  <- (digamma(y+exp(m3)) - digamma(exp(m3)) + log(expit(s3)))*exp(m3)
  dldm3  <- -dldm3
  dlds3  <- (y+exp(m3))*(-expit(s3)) + exp(m3)
  dlds3  <- -dlds3
  
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



if(FALSE){ # Check gradient

 p  <- 3
 L1 <- 5
 L2 <- 2
 n  <- 7
 x  <- matrix(rnorm(n*p),n,p)
 y  <- rbinom(n,10,0.5)

 w   <- init_NegativeBinomial1(p,c(L1,L2))
 test <- predict_NegativeBinomial1(w,x)
 l0  <- loss_NegativeBinomial1(w,x,y)

 numgrad <- function(w,x,y,k,l0,eps){
   g <- w[[k]]
   if(is.vector(g)){
      for(i in 1:length(g)){
        We         <- w
        We[[k]][i] <- We[[k]][i] + eps
        g[i]       <- loss_NegativeBinomial1(We,x,y)
      }
   }
   if(is.matrix(g)){
      for(i in 1:nrow(g)){for(j in 1:ncol(g)){
        We <- w
        We[[k]][i,j] <- We[[k]][i,j] + eps
        g[i,j] <- loss_NegativeBinomial1(We,x,y)
      }}
   }
 return((g-l0)/eps)}

 
 eps <- 10^(-8)
 g1   <- w
 for(k in 1:length(g1)){g1[[k]]<-numgrad(w,x,y,k,l0,eps)}
 g2 <- grad_NegativeBinomial1(w,x,y)
 for(k in 1:length(g1)){
   print(k)
   print(g1[[k]])
   print(g2[[k]])
 }
}


