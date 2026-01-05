
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

predict_hetnorm <- function(w,x){
   m <- act(sweep(x%*%w[[2]],2,w[[1]],"+"))
   m <- act(sweep(m%*%w[[4]],2,w[[3]],"+"))
   m <- m%*%w[[6]] + w[[5]]
   s <- act(sweep(x%*%w[[8]],2,w[[7]],"+"))
   s <- act(sweep(s%*%w[[10]],2,w[[9]],"+"))
   s <- s%*%w[[12]] + w[[11]] 
return(list(mu=m,sigma=exp(s)))}

loss_hetnorm <- function(w,x,y){
   l <- predict_hetnorm(w,x)
   l <- -sum(dnorm(y,l$mu,l$sigma,log=TRUE))
return(l)}

grad_hetnorm <- function(w,x,y){
   M1 <- sweep(x%*%w[[2]],2,w[[1]],"+")
   m1 <- act(M1)
   M2 <- sweep(m1%*%w[[4]],2,w[[3]],"+")
   m2 <- act(M2)
   m3 <- m2%*%w[[6]] + w[[5]]
   S1 <- sweep(x%*%w[[8]],2,w[[7]],"+")
   s1 <- act(S1)
   S2 <- sweep(s1%*%w[[10]],2,w[[9]],"+")
   s2 <- act(S2)
   s3 <- s2%*%w[[12]] + w[[11]] 
   n  <- length(y)

   g  <- w
   
   dldm3  <- -exp(-2*s3)*(y-m3)
   dlds3  <- 1+dldm3*(y-m3)

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

   part2   <- sweep(actp(S2),1,dlds3,"*")
   part1   <- (t(actp(S1))%*%part2)*w[[10]]
   g[[7]]  <- as.vector(part1%*%w[[12]])
   g8      <- NULL
   for(j in 1:ncol(x)){ # remove the loop?
    PP <- sweep(actp(S2),1,dlds3*x[,j],"*")
    PP <- (t(actp(S1))%*%PP)*w[[10]]
    g8 <- rbind(g8,as.vector(PP%*%w[[12]]))
   }
   g[[8]]  <- g8
   g[[9]]  <- colSums(part2)*w[[12]]
   g[[10]] <- sweep(t(s1)%*%part2,2,w[[12]],"*")
   g[[11]] <- colSums(dlds3)
   g[[12]] <- as.vector(t(s2)%*%dlds3)
return(g)}

if(FALSE){ # Check gradient

 p  <- 3
 L1 <- 5
 L2 <- 2
 n  <- 7
 x  <- matrix(rnorm(n*p),n,p)
 y  <- rnorm(n)

 w   <- init_hetnorm(p,L1,L2)
 l0  <- loss_hetnorm(w,x,y)

 numgrad <- function(w,x,y,k,l0,eps){
   g <- w[[k]]
   if(is.vector(g)){
      for(i in 1:length(g)){
        We         <- w
        We[[k]][i] <- We[[k]][i] + eps
        g[i]       <- loss_hetnorm(We,x,y)
      }
   }
   if(is.matrix(g)){
      for(i in 1:nrow(g)){for(j in 1:ncol(g)){
        We <- w
        We[[k]][i,j] <- We[[k]][i,j] + eps
        g[i,j] <- loss_hetnorm(We,x,y)
      }}
   }
 return((g-l0)/eps)}

 
 eps <- 10^(-8)
 g1   <- w
 for(k in 1:12){g1[[k]]<-numgrad(w,x,y,k,l0,eps)}
 g2 <- grad_hetnorm(w,x,y)
 print(g1[[7]])
 print(g2[[7]])
}

