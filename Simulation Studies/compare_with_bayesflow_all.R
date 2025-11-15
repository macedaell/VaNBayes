


# transfer them to the original scale...
dim_vec = dim(simulation_results)
for (i in 1:dim_vec[1]){
  for (j in 1:dim_vec[3]){
    for (k in 1:dim_vec[4]){
      simulation_results[i,,j,k] = original_parameters(simulation_results[i,,j,k])
    }
  }
}



# coverage and MAD
parameter_coverage = matrix(nrow = length(true_parameters), ncol = length(N_vec))
parameter_MAD = matrix(nrow = length(true_parameters), ncol = length(N_vec))
parameter_SDAD = matrix(nrow = length(true_parameters), ncol = length(N_vec))
for (k in 1:length(N_vec)){
  # coverage
  parameter_coverage[,k] = c(sum((simulation_results[,1,1,k] < true_parameters[1]) & (true_parameters[1] < simulation_results[,1,3,k]))/nsims,
                             sum((simulation_results[,2,1,k] < true_parameters[2]) & (true_parameters[2] < simulation_results[,2,3,k]))/nsims,
                             sum((simulation_results[,3,1,k] < true_parameters[3]) & (true_parameters[3] < simulation_results[,3,3,k]))/nsims)
  # MAD
  parameter_MAD[,k] = c(mean(abs(simulation_results[,1,2,k] - true_parameters[1])), 
                        mean(abs(simulation_results[,2,2,k] - true_parameters[2])), 
                        mean(abs(simulation_results[,3,2,k] - true_parameters[3])))
  
  # SD
  parameter_SDAD[,k] = c(sd(abs(simulation_results[,1,2,k] - true_parameters[1])), 
                         sd(abs(simulation_results[,2,2,k] - true_parameters[2])), 
                         sd(abs(simulation_results[,3,2,k] - true_parameters[3])))
  
  
}

print("VaNBayes Parameter Coverage:")
print(parameter_coverage)
print("VaNBayes Parameter MAD")
print(parameter_MAD)
print("VaNBayes Parameter SDAD")
print(parameter_SDAD)

VaNBayes_medians = simulation_results[,,2,]


# import Bayesflow
bayesflow_medians = array(0, dim = dim(VaNBayes_medians))
for (k in 1:length(N_vec)){
  N = format(N_vec[k], scientific = FALSE)
  bayesflow_medians[,,k] = as.matrix(read.csv(paste0("./bayesflow_posterior_medians_N_",N,".csv"), header = FALSE))
}
bayesflow_parameter_coverage = read.csv("./bayesflow_parameter_coverage.csv", header = FALSE)
bayesflow_parameter_coverage = as.matrix(bayesflow_parameter_coverage[,1:length(N_vec)])
bayesflow_MAD = read.csv("./bayesflow_parameter_MAD.csv", header = FALSE)
bayesflow_MAD = as.matrix(bayesflow_MAD[,1:length(N_vec)])
bayesflow_SDAD = read.csv("./bayesflow_parameter_SDAD.csv", header = FALSE)
bayesflow_SDAD = as.matrix(bayesflow_SDAD[,1:length(N_vec)])





## MAD TRAJECTORY ####
for (j in 1:length(true_parameters)){ # true lambda, I0, psi
  plot(x = 1:length(N_vec),
       y = parameter_MAD[j,],
       type = "l",
       lwd = 2, 
       lty = 1,
       axes = F,
       cex.lab=1.5,
       main = parameter_symbols[j],
       col = "lightblue",
       xlab = "Number of Simulated Samples",
       ylab = "MAD",
       xlim = range(c(1:length(N_vec)-0.25, 1:length(N_vec)+0.25)),
       ylim = range(c(parameter_MAD[j,] - 2*parameter_SDAD[j,], 
                      parameter_MAD[j,] + 2*parameter_SDAD[j,], 
                      bayesflow_MAD[j,] - 2*bayesflow_SDAD[j,],
                      bayesflow_MAD[j,] + 2*bayesflow_SDAD[j,])))
  lines(x = 1:length(N_vec),
        y = parameter_MAD[j,] - 2*parameter_SDAD[j,],
        lwd = 2, 
        lty = 2,
        col = "lightblue",)
  lines(x = 1:length(N_vec),
        y = parameter_MAD[j,] + 2*parameter_SDAD[j,],
        lwd = 2, 
        lty = 2,
        col = "lightblue",)
  lines(x = 1:length(N_vec),
        y = bayesflow_MAD[j,],
        lwd = 2, 
        lty = 1,
        col = "darkblue",)
  lines(x = 1:length(N_vec),
        y = bayesflow_MAD[j,] - 2*bayesflow_SDAD[j,],
        lwd = 2, 
        lty = 2,
        col = "darkblue",)
  lines(x = 1:length(N_vec),
        y = bayesflow_MAD[j,] + 2*bayesflow_SDAD[j,],
        lwd = 2, 
        lty = 2,
        col = "darkblue",)
  abline(h = 0, col = "red")
  axis(2,cex.axis=1.25,cex.lab=1.25)
  axis(1,at=1:length(N_vec),lab=N_vec*(1-PROPVAL),cex.axis=1.5,cex.lab=1.5)
  legend("topright", inset=c(0,-.1),xpd=TRUE, # need to adjust the first inset number
         legend = c("VaNBayes", "Bayesflow"),
         fill=c("lightblue", "darkblue"),
         bty="n",cex=1.5)
}




## Boxplots ####
# Boxplot with scaled x-axis
for (j in 1:length(true_parameters)){ # true lambda, I0, psi
  spacing = 0.1
  plot(x = 1:length(N_vec),
       y = rep(NA,length(N_vec)),
       axes = F,
       cex.lab=1.5,
       main = parameter_symbols[j],
       xlab = "Number of Simulated Samples",
       ylab = "Posterior Medians",
       xlim = range(c(1:length(N_vec)-0.25, 1:length(N_vec)+0.25)),
       # xlim = range(c(min(N_vec)-0.25, max(N_vec)+0.25)),
       ylim = range(c(VaNBayes_medians[,j,], bayesflow_medians[,j,])))
  abline(h = 0, col = "black", lwd = 2)
  for (k in 1:length(N_vec)){
    boxplot(VaNBayes_medians[,j,k],at=k-spacing,add=TRUE,axes=F,col="lightblue",boxwex=.4,outline=FALSE)
    boxplot(bayesflow_medians[,j,k],at=k+spacing,add=TRUE,axes=F,col="darkblue",boxwex=.4,outline=FALSE)
  }
  
  # boxplot(VaNBayes_medians[,j,1],at=1-spacing,add=TRUE,axes=F,col="lightblue",boxwex=.4,outline=FALSE)
  # boxplot(bayesflow_medians[,j,1],at=1+spacing,add=TRUE,axes=F,col="darkblue",boxwex=.4,outline=FALSE)
  # boxplot(VaNBayes_medians[,j,2],at=2-spacing,add=TRUE,axes=F,col="lightblue",boxwex=.4,outline=FALSE)
  # boxplot(bayesflow_medians[,j,2],at=2+spacing,add=TRUE,axes=F,col="darkblue",boxwex=.4,outline=FALSE)
  # boxplot(VaNBayes_medians[,j,3],at=3-spacing,add=TRUE,axes=F,col="lightblue",boxwex=.4,outline=FALSE)
  # boxplot(bayesflow_medians[,j,3],at=3+spacing,add=TRUE,axes=F,col="darkblue",boxwex=.4,outline=FALSE)
  # boxplot(VaNBayes_medians[,j,4],at=4-spacing,add=TRUE,axes=F,col="lightblue",boxwex=.4,outline=FALSE)
  # boxplot(bayesflow_medians[,j,4],at=4+spacing,add=TRUE,axes=F,col="darkblue",boxwex=.4,outline=FALSE)
  # boxplot(VaNBayes_medians[,j,5],at=5-spacing,add=TRUE,axes=F,col="lightblue",boxwex=.4,outline=FALSE)
  # boxplot(bayesflow_medians[,j,5],at=5+spacing,add=TRUE,axes=F,col="darkblue",boxwex=.4,outline=FALSE)
  # boxplot(VaNBayes_medians[,j,6],at=6-spacing,add=TRUE,axes=F,col="lightblue",boxwex=.4,outline=FALSE)
  # boxplot(bayesflow_medians[,j,6],at=6+spacing,add=TRUE,axes=F,col="darkblue",boxwex=.4,outline=FALSE)
  # boxplot(VaNBayes_medians[,j,7],at=7-spacing,add=TRUE,axes=F,col="lightblue",boxwex=.4,outline=FALSE)
  # boxplot(bayesflow_medians[,j,7],at=7+spacing,add=TRUE,axes=F,col="darkblue",boxwex=.4,outline=FALSE)
  abline(h = true_parameters[j], col = "red")
  axis(2,cex.axis=1.25,cex.lab=1.25)
  axis(1,at=1:length(N_vec),lab=N_vec*(1-PROPVAL),cex.axis=1.5,cex.lab=1.5)
  # axis(1,at=N_vec,lab=N_vec,cex.axis=1.5,cex.lab=1.5)
  legend("topright", inset=c(0,-.1),xpd=TRUE, # need to adjust the first inset number
         legend = c("VaNBayes", "Bayesflow"),
         fill=c("lightblue", "darkblue"),
         bty="n",cex=1.5)
}

