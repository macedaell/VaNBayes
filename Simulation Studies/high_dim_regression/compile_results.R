



# cases of p
p_vec = c(5, 10, 15 , 20, 25, 30, 35, 40)



# line graph for MAE's -- include average betas, sigma
VaNBayes_beta_MAE = numeric(length(p_vec))
VaNBayes_sigma_MAE = numeric(length(p_vec))
Bayesflow_beta_MAE = numeric(length(p_vec))
Bayesflow_sigma_MAE = numeric(length(p_vec))

# uncertainty quantification
VaNBayes_beta_MAE_lower = numeric(length(p_vec))
VaNBayes_beta_MAE_upper = numeric(length(p_vec))
Bayesflow_beta_MAE_lower = numeric(length(p_vec))
Bayesflow_beta_MAE_upper = numeric(length(p_vec))

# line graph for coverages -- average over all betas, sigma
VaNBayes_beta_COV = numeric(length(p_vec))
VaNBayes_sigma_COV = numeric(length(p_vec))
Bayesflow_beta_COV = numeric(length(p_vec))
Bayesflow_sigma_COV = numeric(length(p_vec))

# computing all the values
for (i in 1:length(p_vec)){
  
  p = p_vec[i]
  
  # beta values
  beta0 = 2*((1:p %% 5) - 2)
  
  
  print(paste0("Loading for p=", p))
  
  # import VaNBayes
  load(paste0("./results/VaNBayes_results_p", p, ".RData"))
  
  # obtain VaNBayes coverages
  VaNBayes_sigma_coverage = sigma_coverage
  VaNBayes_alpha_coverage = alpha_coverage
  VaNBayes_beta_coverage = beta_coverage
  
  # obtain VaNBayes MADs
  VaNBayes_sigma_MAD = sigma_MAD
  VaNBayes_alpha_MAD = alpha_MAD
  VaNBayes_beta_MAD = beta_MAD
  
  # obtain VaNBayes medians
  VaNBayes_sigma_medians = sigma_medians
  VaNBayes_alpha_medians = alpha_medians
  VaNBayes_beta_medians = beta_medians
  
  VaNBayes_individual_MAEs = abs(sweep(VaNBayes_beta_medians, 2, beta0, "-"))
  VaNBayes_this_p_MAEs = apply(VaNBayes_individual_MAEs, 1, mean)
  
  
  # import Bayesflow
  bayesflow_coverage = read.csv(paste0("./high_dim_regression_bayesflow/results/bayesflow_parameter_coverage_p",p,".csv"), header = FALSE)
  bayesflow_MAD = read.csv(paste0("./high_dim_regression_bayesflow/results/bayesflow_parameter_MAD_p",p,".csv"), header = FALSE)
  bayesflow_medians = read.csv(paste0("./high_dim_regression_bayesflow/results/bayesflow_posterior_medians_p",p,".csv"), header = FALSE)
  
  # obtain Bayesflow coverages
  Bayesflow_sigma_coverage = bayesflow_coverage[1,]
  Bayesflow_alpha_coverage = bayesflow_coverage[2,]
  Bayesflow_beta_coverage = bayesflow_coverage[-c(1,2),]
  
  # obtain Bayesflow MADs
  Bayesflow_sigma_MAD = bayesflow_MAD[1,]
  Bayesflow_alpha_MAD = bayesflow_MAD[2,]
  Bayesflow_beta_MAD = bayesflow_MAD[-c(1,2),]
  
  # obtain Bayesflow medians
  Bayesflow_sigma_medians = bayesflow_medians[,1]
  Bayesflow_alpha_medians = bayesflow_medians[,2]
  Bayesflow_beta_medians = bayesflow_medians[,-c(1,2)]
  
  Bayesflow_individual_MAEs = abs(sweep(Bayesflow_beta_medians, 2, beta0, "-"))
  Bayesflow_this_p_MAEs = apply(Bayesflow_individual_MAEs, 1, mean)
  
  
  # new MAEs-- also won't use the Sigma's
  VaNBayes_beta_MAE[i] = mean(VaNBayes_this_p_MAEs)
  VaNBayes_beta_MAE_se = sd(VaNBayes_this_p_MAEs)/sqrt(100)
  VaNBayes_beta_MAE_lower[i] = VaNBayes_beta_MAE[i] - 2*VaNBayes_beta_MAE_se
  VaNBayes_beta_MAE_upper[i] = VaNBayes_beta_MAE[i] + 2*VaNBayes_beta_MAE_se
  Bayesflow_beta_MAE[i] = mean(Bayesflow_this_p_MAEs)
  Bayesflow_beta_MAE_se = sd(VaNBayes_this_p_MAEs)/sqrt(100)
  Bayesflow_beta_MAE_lower[i] = Bayesflow_beta_MAE[i] - 2*Bayesflow_beta_MAE_se
  Bayesflow_beta_MAE_upper[i] = Bayesflow_beta_MAE[i] + 2*Bayesflow_beta_MAE_se
  
  # we won't use these
  VaNBayes_sigma_MAE[i] = mean(VaNBayes_sigma_MAD)
  Bayesflow_sigma_MAE[i] = mean(Bayesflow_sigma_MAD)
  
  
  # # old MAEs
  # VaNBayes_beta_MAE[i] = mean(VaNBayes_beta_MAD)
  # VaNBayes_sigma_MAE[i] = mean(VaNBayes_sigma_MAD)
  # Bayesflow_beta_MAE[i] = mean(Bayesflow_beta_MAD)
  # Bayesflow_sigma_MAE[i] = mean(Bayesflow_sigma_MAD)
  
  # old
  VaNBayes_beta_COV[i] = mean(VaNBayes_beta_coverage)
  VaNBayes_sigma_COV[i] = mean(VaNBayes_sigma_coverage)
  Bayesflow_beta_COV[i] = mean(Bayesflow_beta_coverage)
  Bayesflow_sigma_COV[i] = mean(Bayesflow_sigma_coverage)
  
  # line plot of betas, for each p
  all_plot_values = c(VaNBayes_beta_MAD, Bayesflow_beta_MAD)
  plot(x = 1:p, y = VaNBayes_beta_MAD, type = "l", col = "black", lty = 1, lwd = 2, xlab = "number of covariates", ylab = "MAE", ylim = c(min(all_plot_values), max(all_plot_values)))
  lines(x = 1:p,, y = Bayesflow_beta_MAD, col = "red", lty = 1, lwd = 2)
  legend("topleft",c("VaNBayes", "BayesFlow"),fill=c("black","red"),bty="n",cex=1.5)
  
  
  cat("Bayesflow first: ", Bayesflow_beta_MAD[1], "--Bayesflow last: ", Bayesflow_beta_MAD[p], "\n")
  cat("VaNBayes first: ", VaNBayes_beta_MAD[1], "--VaNBayes last: ", VaNBayes_beta_MAD[p], "\n")
  
  
  # combined boxplot
  plot(NA,xlim=c(0,p+1), ylim=range(c(Bayesflow_beta_medians, VaNBayes_beta_medians)),axes=F,ylab="Posterior Median", xlab = "Covariate",cex.lab=1.5)
  boxplot(VaNBayes_beta_medians,at=1:p-0.25,add=TRUE,axes=F,col="white",boxwex=.4,outline=FALSE)
  boxplot(Bayesflow_beta_medians,at=1:p+0.25,add=TRUE,axes=F,col= "darkgray",boxwex=.4,outline=FALSE)
  points(x = 1:p, y = beta0, col=2,pch=19)
  axis(2,cex.axis=1.25,cex.lab=1.25)
  axis(1,at=1:p,lab=1:p,cex.axis=1.25,cex.lab=1.25)
  legend("topright",c("VaNBayes", "BayesFlow"),fill=c("white","darkgray"),bty="n",cex=1.4)
  
  
}

# line plot of betas
all_plot_values = c(VaNBayes_beta_MAE, Bayesflow_beta_MAE)
plot(x = p_vec, y = VaNBayes_beta_MAE, type = "l", col = "black", lty = 1, lwd = 2, xlab = "Number Of Covariates", ylab = "Mean Absolute Error", ylim = c(min(all_plot_values), max(all_plot_values)),cex.lab=1.5, cex.axis = 1.5)
lines(x = p_vec, y = VaNBayes_beta_MAE_lower, col = "black", lty = 2, lwd = 2)
lines(x = p_vec, y = VaNBayes_beta_MAE_upper, col = "black", lty = 2, lwd = 2)
lines(x = p_vec, y = Bayesflow_beta_MAE, col = "red", lty = 1, lwd = 2)
lines(x = p_vec, y = Bayesflow_beta_MAE_lower, col = "red", lty = 2, lwd = 2)
lines(x = p_vec, y = Bayesflow_beta_MAE_upper, col = "red", lty = 2, lwd = 2)
legend("topleft",c("VaNBayes", "BayesFlow"),fill=c("black","red"),bty="n",cex=1.5)

# line plot of sigmas
all_plot_values = c(VaNBayes_sigma_MAE, Bayesflow_sigma_MAE)
plot(x = p_vec, y = VaNBayes_sigma_MAE, type = "l", col = "black", lty = 1, lwd = 2, xlab = "Number Of Covariates", ylab = "MAE", ylim = c(min(all_plot_values), max(all_plot_values)), main = "SIGMA")
lines(x = p_vec, y = Bayesflow_sigma_MAE, col = "red", lty = 1, lwd = 2)
legend("topleft",c("VaNBayes", "BayesFlow"),fill=c("black","red"),bty="n",cex=1.5)

# line plot of beta-coverages
all_plot_values = c(VaNBayes_beta_COV, Bayesflow_beta_COV, .90)
plot(x = p_vec, y = VaNBayes_beta_COV, type = "l", col = "black", lty = 1, lwd = 2, xlab = "Number Of Covariates", ylab = "COV", ylim = c(min(all_plot_values), max(all_plot_values)), main = "BETAS")
lines(x = p_vec, y = Bayesflow_beta_COV, col = "red", lty = 1, lwd = 2)
abline(h = 0.90, col = "blue", lwd = 2)
legend("topleft",c("VaNBayes", "BayesFlow"),fill=c("black","red"),bty="n",cex=1.5)

# line plot of sigma-coverages
all_plot_values = c(VaNBayes_sigma_COV, Bayesflow_sigma_COV, .90)
plot(x = p_vec, y = VaNBayes_sigma_COV, type = "l", col = "black", lty = 1, lwd = 2, xlab = "Number Of Covariates", ylab = "COV", ylim = c(min(all_plot_values), max(all_plot_values)), main = "SIGMA")
lines(x = p_vec, y = Bayesflow_sigma_COV, col = "red", lty = 1, lwd = 2)
abline(h = 0.90, col = "blue", lwd = 2)
legend("topleft",c("VaNBayes", "BayesFlow"),fill=c("black","red"),bty="n",cex=1.5)



