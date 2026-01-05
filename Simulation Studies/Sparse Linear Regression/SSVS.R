SSVS <- function(y,x,xp=NULL,init_beta=NULL,init_int=NULL,init_sigma2=NULL,
                 iters=20000,burn=5000,update=Inf,
                 a=0.1,b=0.1,c=100000,d=100000,tau=1){

    n  <- length(y)
    p  <- ncol(x)
    np <- ifelse(is.null(xp),0,nrow(xp))

    #initial values
    int    <- mean(y)
    beta   <- alpha <- inout<-rep(0,p)
    sigma2 <- var(y)
    p_in   <- c/(c+d)
    if(!is.null(init_beta)){
        alpha  <- init_beta
        inout  <- ifelse(beta==0,0,1)
        beta   <- inout*alpha
    }
    if(!is.null(init_int)){
        int <- init_int
    }
    if(!is.null(init_sigma2)){
        sigma2 <- init_sigma2
    }


    #keep track of stuff:
    keep.beta           <- matrix(0,iters,p)
    keep.sigma2         <- rep(0,iters)
    keep.int            <- rep(0,iters)
    keep.p_in           <- rep(0,iters)
    keep.pred           <- matrix(0,iters,np)

    for(i in 1:iters){

      sigma2 <- 1/rgamma(1,n/2+a,sum((y-int-x%*%beta)^2)/2+b)

      VVV <- n/sigma2 + 1/tau^2
      MMM <- sum(y-x%*%beta)/sigma2
      int <- rnorm(1,MMM/VVV,1/sqrt(VVV))
 
      #update alpha
      z     <- x%*%diag(inout)
      V     <- solve(t(z)%*%z/sigma2+diag(p)/tau^2)
      M     <- t(z)%*%(y-int)/sigma2
      alpha <- V%*%M+t(chol(V))%*%rnorm(p)
      beta  <- alpha*inout

      #update inclusion indicators: 
      r<-y-int-x%*%beta
      for(j in 1:p){
        r         <- r+x[,j]*beta[j]
        log.p.in  <- log(p_in)-0.5*sum((r-x[,j]*alpha[j])^2)/sigma2
        log.p.out <- log(1-p_in)-0.5*sum(r^2)/sigma2
        diff      <- log.p.in-log.p.out
        if(diff>10){diff<-10}
        p1        <- exp(diff)/(1+exp(diff))
        inout[j]  <- rbinom(1,1,p1)
        beta[j]   <- inout[j]*alpha[j]
        r         <- r-x[,j]*beta[j]
      }

      p_in <- rbeta(1,sum(inout==1)+c,sum(inout==0)+d)

      # Keep track of stuff:
      keep.beta[i,]  <- beta
      keep.sigma2[i] <- sigma2
      keep.int[i]    <- int
      keep.p_in[i]   <- p_in
      if(np>0){
        keep.pred[i,] <- rnorm(np,int + xp%*%beta,sqrt(sigma2))
      }

      if(i%%update==0){plot(beta,main=paste("Iteration",i))}
    }
    output <- list(beta=keep.beta,sigma=sqrt(keep.sigma2),
                   int=keep.int,p_in=keep.p_in,
                   pred=keep.pred)
return(output)}