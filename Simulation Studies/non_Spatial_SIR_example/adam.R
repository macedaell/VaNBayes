adam <- function(w, x, y, loss, grad,
                 lr = 0.001, batchsize = 100, epochs = 50, propval = 0.2,
                 beta1=0.9, beta2=0.999, epsilon=0.1^8,
                 verbose=2,early_stop=5){
  
  # https://arxiv.org/pdf/1412.6980.pdf
  
  val   <- seq(0,1,length=length(y)) > 1-propval
  x_val <- x[val,]
  y_val <- y[val]
  x     <- x[!val,]
  y     <- y[!val]
  n     <- length(y)
  nb    <- floor(n/batchsize)
  bs    <- rep(1:nb,length=n)[1:n]
  
  m <- v <- w
  for(k in 1:length(w)){
    m[[k]] <- 0*m[[k]]
    v[[k]] <- 0*v[[k]]
  }  
  
  train_loss <- val_loss <- rep(NA,epochs)
  
  t <- 1
  stop = FALSE
  for(iter in 1:epochs){if(!stop){
    batch <- sample(bs)
    w0  <- w; m0  <- m; v0  <- v
    for(b in 1:nb){
      sub <- which(batch==b)
      N   <- n/length(sub)
      g   <- grad(w,x[sub,],y[sub])
      for(k in 1:length(w)){
        W         <- w        
        G         <- N*g[[k]]
        m[[k]]    <- beta1*m[[k]] + (1-beta1)*G
        v[[k]]    <- beta2*v[[k]] + (1-beta2)*(G^2)
        mhat      <- m[[k]]/(1-(beta1^t))
        vhat      <- v[[k]]/(1-(beta2^t))
        w[[k]]    <- w[[k]] - lr*mhat/(sqrt(vhat)+epsilon) 
        if(sum(is.na(w[[k]]))>0){w<-W;lr<-lr/2}
      } 
      t <- t+1
    }
    if(is.na(loss(w,x,y))){w<-w0;m<-m0;v<-v0}
    if(is.na(loss(w,x_val,y_val))){w<-w0;m<-m0;v<-v0}
    
    train_loss[iter] <- loss(w,x,y)/length(y)
    val_loss[iter]   <- loss(w,x_val,y_val)/length(y_val)
    if(val_loss[iter]==min(val_loss[1:iter])){bestw<-w}   
    
    if(verbose == 1){
      print(c(b,train_loss[iter],val_loss[iter]))
    }
    if(verbose == 2){
      plot(train_loss,xlab="Epoch",type="l",ylab="Loss",col=1,xlim=c(1,iter),
           ylim=range(c(val_loss,train_loss),na.rm=TRUE))
      lines(val_loss,col=2)
      legend("topright",c("Training","Validation"),lty=1,col=1:2,bty="n")     
    }
    
    if(iter>early_stop){
      j    <- iter-early_stop
      stop <- min(val_loss[which(1:iter > j)]-val_loss[j])>0        
    }
  }}   
  tuning <- list(lr = lr, batchsize = batchsize, epochs = epochs, propval = propval,
                 beta1=beta1, beta2=beta2, epsilon=epsilon, early_stop=early_stop)
  output <- list(w=bestw,tuning=tuning,train_loss=train_loss,val_loss=val_loss)
  return(output)}

if(FALSE){
  n <- 1000
  p <- 10
  x <- matrix(rnorm(n*p),n,p)
  y <- rnorm(n,x[,1]-x[,2],1)
  print("True optimum")
  print(round(lm(y~x)$coef,3))
  
  loss_ols <- function(w,x,y){
    beta <- w[[1]]
    X    <- cbind(1,x)
    sse  <- sum((y-X%*%beta)^2)
    return(sse)}
  
  grad_ols <- function(w,x,y){
    beta <- w[[1]]
    X    <- cbind(1,x)
    g    <- list(-2*as.vector(t(X)%*%(y-X%*%beta)))
    return(g)}
  
  
  init <- list(rnorm(p+1))
  print("Initial value")
  print(round(init[[1]],3))
  fit  <- adam(w=init, x=x, y=y, loss=loss_ols, grad=grad_ols,
               lr=0.001,epochs=100)
  print("Estimate")
  print(round(fit$w[[1]],3))
}

# Functions to scale the inputs:

stdq <- function(x,q,low=-1,high=1){
  m         <- length(q)
  x[x<q[1]] <- q[1]
  x[x>q[m]] <- q[m]
  U         <- 1+0*x
  for(j in 2:(m-1)){U[x>q[j]]<-j} 
  U         <- U+(x-q[U])/(q[U+1]-q[U])
  return(low + (high-low)*U/m)}

quantile_scale <- function(x,q){
  for(k in 1:ncol(x)){
    x[,k] <- stdq(x[,k],q[,k])
  }
  return(x)}
