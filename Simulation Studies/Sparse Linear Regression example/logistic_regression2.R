
logit <- function(x){log(x)-log(1-x)}
expit <- function(x){1/(1+exp(-x))}
act   <- function(x){pmax(x,0)}
actp  <- function(x){x>0}

init_logistic <- function(p,L1,L2,init_mn=0,sigma=0.01){
 w       <- list()
 w[[1]]  <- sigma*rnorm(L1)
 w[[2]]  <- sigma*matrix(rnorm(p*L1),p,L1)
 w[[3]]  <- sigma*rnorm(L2)
 w[[4]]  <- sigma*matrix(rnorm(L1*L2),L1,L2)
 w[[5]]  <- init_mn + sigma*rnorm(1)
 w[[6]]  <- sigma*rnorm(L2)
return(w)}

predict_logistic <- function(w,x){
   m <- act(sweep(x%*%w[[2]],2,w[[1]],"+"))
   m <- act(sweep(m%*%w[[4]],2,w[[3]],"+"))
   m <- expit(m%*%w[[6]] + w[[5]])
return(m)}

loss_logistic <- function(w,x,y){
   l <- predict_logistic(w,x)
   l <- -sum(y*log(l)+(1-y)*log(1-l))
return(l)}

grad_logistic <- function(w,x,y){
   M1 <- sweep(x%*%w[[2]],2,w[[1]],"+")
   m1 <- act(M1)
   M2 <- sweep(m1%*%w[[4]],2,w[[3]],"+")
   m2 <- act(M2)
   m3 <- m2%*%w[[6]] + w[[5]]
   n  <- length(y)

   g  <- w
   
   dldm3  <- expit(m3)-y
   part2  <- sweep(actp(M2),1,dldm3,"*")
   part1  <- (t(actp(M1))%*%part2)*w[[4]]
   g[[1]] <- as.vector(part1%*%w[[6]])
   g2     <- NULL
   for(j in 1:ncol(x)){ # remove the loop?
    PP <- sweep(actp(M2),1,dldm3*x[,j],"*")
    PP <- (t(actp(M1))%*%PP)*w[[4]]
    g2 <- rbind(g2,as.vector(PP%*%w[[6]]))
   }
   g[[2]] <- g2
   g[[3]] <- colSums(part2)*w[[6]]
   g[[4]] <- sweep(t(m1)%*%part2,2,w[[6]],"*")
   g[[5]] <- colSums(dldm3)
   g[[6]] <- as.vector(t(m2)%*%dldm3)
return(g)}


if(FALSE){ # Check gradient

 p  <- 3
 L1 <- 5
 L2 <- 2
 n  <- 7
 x  <- matrix(rnorm(n*p),n,p)
 y  <- rbinom(n,1,pnorm(x[,3]))

 w   <- init_logistic(p,L1,L2)
 l0  <- loss_logistic(w,x,y)

 numgrad2 <- function(w,x,y,k,l0,eps){
   g <- w[[k]]
   if(is.vector(g)){
      for(i in 1:length(g)){
        We         <- w
        We[[k]][i] <- We[[k]][i] + eps
        g[i]       <- loss_logistic(We,x,y)
      }
   }
   if(is.matrix(g)){
      for(i in 1:nrow(g)){for(j in 1:ncol(g)){
        We <- w
        We[[k]][i,j] <- We[[k]][i,j] + eps
        g[i,j] <- loss_logistic(We,x,y)
      }}
   }
 return((g-l0)/eps)}

 
 eps <- 10^(-8)
 g1   <- w
 for(k in 1:6){g1[[k]]<-numgrad2(w,x,y,k,l0,eps)}
 g2 <- grad_logistic(w,x,y)
 for(j in 1:6){
   print(paste0("w",j))
   print(g1[[j]])
   print(g2[[j]])
 }
}
