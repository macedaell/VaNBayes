library(fields)
library(spam)
logit <- function(x){log(x)-log(1-x)}
expit <- function(x){1/(1+exp(-x))}

###############################################################
# Draw samples from the distribution in Equation (7) of
# https://link.springer.com/article/10.1007/s13253-014-0183-0
###############################################################

# Here N is a list
rautologistic <- function(kappa,gamma,N,w=NULL,iters=50){
   n  <- length(N)
   lk <- logit(kappa)  
   if(is.null(w)){w<-rbinom(n,1,kappa)}
   lp <- lk
   for(i in 1:n){
        lp[i] <- lk[i] - gamma*sum(kappa[N[[i]]])
   }  
   for(iter in 1:iters){for(i in 1:n){
     w[i] <- logit(runif(1)) < lp[i] + gamma*sum(w[N[[i]]])
   }}
return(w)}

# Here N is a sparse matrix
rautologisticblock <- function(kappa,gamma,N,block=NULL,w=NULL,iters=100){
   n   <- nrow(N)
   if(is.null(w)){w<-rbinom(n,1,kappa)}
   lp  <- logit(kappa) - gamma*N%*%kappa
   b1  <- which(block==1)
   b2  <- which(block==2)
   lp1 <- lp[b1]
   lp2 <- lp[b2]
   N1  <- N[b1,]
   N2  <- N[b2,]
   n1  <- length(b1)
   n2  <- length(b2)
   for(iter in 1:iters){
     w[b1] <- logit(runif(n1)) < lp1 + gamma*N1%*%w
     w[b2] <- logit(runif(n2)) < lp2 + gamma*N2%*%w
   }
return(w)}

# Function for feature construction
moranI <- function(y,N){
   den <- var(y)
   num <- 0
   for(i in 1:length(y)){
      num <- num + y[i]*sum(y[N[[i]]])
   }
   num <- num/sum(unlist(lapply(N,length)))
return(num/den)}

GearyC <- function(y,N){
   den <- var(y)
   num <- 0
   for(i in 1:length(y)){
      num <- num + sum((y[i]-y[N[[i]]])^2)
   }
   num <- num/sum(unlist(lapply(N,length)))
return(0.5*num/den)}


