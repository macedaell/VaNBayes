rm(list=ls())
setwd("S:\\Desktop\\SparseLinModel\\Sim\\")

############################################
##########          PIPs          ##########
############################################

pip_box <- function(p0,p1){
  p <- ncol(p0)
  plot(NA,xlim=c(0.5,p+0.5),ylim=0:1,axes=F,xlab="Covariate",ylab="Inclusion probability",
       cex.lab=1.5)
  for(j in 1:p){
    boxplot(p0[,j],at=j-0.15,add=TRUE,axes=F,col="white",boxwex=.5,outline=FALSE)
    boxplot(p1[,j],at=j+0.15,add=TRUE,axes=F,col= "gray",boxwex=.5,outline=FALSE)
  }
  axis(1,at=1:p,cex.axis=1.5)
  axis(2,at=seq(0,1,.2),cex.axis=1.5)
  legend("topright",c("MCMC","VaNBayes"),fill=c("white","gray"),bty="n",cex=1.5)
}

output <- NULL
for(p in c(10,20)){for(type in 1:6){
  load(paste0("Sim_p",p,"fit",type,".RData"))
  out <- c(p,nodes,mean(CE),mean(CA),mean(BS),cor(as.vector(pip0),as.vector(pip1)))
  output <- rbind(output,out)
  pdf(paste0("Sim_p",p,"fit",type,".pdf"))
  pip_box(pip0,pip1)
  dev.off()
}}

out <- cbind(output[1:6,2:6],output[7:12,4:6])
round(out,4)

###########################################
############       Sigma       ############
###########################################

for(p in c(10,20)){for(jjj in 1:6){
  load(paste0("Sim_p",p,"fit",jjj,".RData"))
  print(c(p,jjj,LS_sig))
}}

type  <- 4
p     <- 10
load(paste0("Sim_p",p,"fit",type,".RData"))
cov0  <- mean((q_sig0[,1]<sigma0) & (sigma0<q_sig0[,3]))
cov1  <- mean((q_sig1[,1]<sigma0) & (sigma0<q_sig1[,3]))
mad0  <- mean(abs(q_sig0[,2]-sigma0))
mad1  <- mean(abs(q_sig1[,2]-sigma0))
out   <- matrix(c(mad0,cov0,mad1,cov1),2,2,byrow=TRUE)
colnames(out) <- c("MSE","COV")
rownames(out) <- c("MCMC","VaNBayes")
print(round(out,3))

s0_10 <-q_sig0[,2]
s1_10 <-q_sig1[,2]
 

type  <- 2
p     <- 20
load(paste0("Sim_p",p,"fit",type,".RData"))
cov0  <- mean((q_sig0[,1]<sigma0) & (sigma0<q_sig0[,3]))
cov1  <- mean((q_sig1[,1]<sigma0) & (sigma0<q_sig1[,3]))
mad0  <- mean(abs(q_sig0[,2]-sigma0))
mad1  <- mean(abs(q_sig1[,2]-sigma0))
out   <- matrix(c(mad0,cov0,mad1,cov1),2,2,byrow=TRUE)
colnames(out) <- c("MSE","COV")
rownames(out) <- c("MCMC","VaNBayes")
print(round(out,3))

s0_20 <-q_sig0[,2]
s1_20 <-q_sig1[,2]

pdf("sigma_medians.pdf")
 r <- range(c(s0_10,s1_10,s0_20,s1_20))
 plot(s0_10,s1_10,xlab="MCMC",ylab="VaNBayes",xlim=r,ylim=r)
 abline(0,1)
 plot(s0_20,s1_20,xlab="MCMC",ylab="VaNBayes",xlim=r,ylim=r)
 abline(0,1)
dev.off()


###########################################
############    Predictions    ############
###########################################

for(p in c(10,20)){for(jjj in 1:6){
  load(paste0("Sim_p",p,"fit",jjj,".RData"))
  print(c(p,jjj,LS_pred))
}}


type  <- 5
p     <- 10
load(paste0("Sim_p",p,"fit",type,".RData"))
cov0  <- cov1 <- mad0 <- mad1 <- 0
for(j in 1:np){
  cov0  <- cov0 + mean((q_pred0[,j,1]<Y_p[,j]) & (Y_p[,j]<q_pred0[,j,3]))/np
  cov1  <- cov1 + mean((q_pred1[,j,1]<Y_p[,j]) & (Y_p[,j]<q_pred1[,j,3]))/np
  mad0  <- mad0 + mean(abs(q_pred0[,j,2]-Y_p[,j]))/np
  mad1  <- mad1 + mean(abs(q_pred1[,j,2]-Y_p[,j]))/np
}
out   <- matrix(c(mad0,cov0,mad1,cov1),2,2,byrow=TRUE)
colnames(out) <- c("MSE","COV")
rownames(out) <- c("MCMC","VaNBayes")
print(round(out,3))


type  <- 6
p     <- 20
load(paste0("Sim_p",p,"fit",type,".RData"))
cov0  <- cov1 <- mad0 <- mad1 <- 0
for(j in 1:np){
  cov0  <- cov0 + mean((q_pred0[,j,1]<Y_p[,j]) & (Y_p[,j]<q_pred0[,j,3]))/np
  cov1  <- cov1 + mean((q_pred1[,j,1]<Y_p[,j]) & (Y_p[,j]<q_pred1[,j,3]))/np
  mad0  <- mad0 + mean(abs(q_pred0[,j,2]-Y_p[,j]))/np
  mad1  <- mad1 + mean(abs(q_pred1[,j,2]-Y_p[,j]))/np
}
out   <- matrix(c(mad0,cov0,mad1,cov1),2,2,byrow=TRUE)
colnames(out) <- c("MSE","COV")
rownames(out) <- c("MCMC","VanBayes")
print(round(out,3))


