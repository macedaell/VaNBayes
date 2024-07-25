adam.R: The optimization algorithm for the neural network.
hetnorm2.R: Includes the initialization, loss function, and predictions for our variational posterior related to the sigma parameter and prediction.
logistic_regression2.R: Includes the initialization, loss function, and predictions for our posterior related to the logistic regression weights.
SSVS.R: Includes the function that performs MCMC to compute the PIP's and posterior estimates of the parameters.
Sim.R: The source file. First, generates datasets from our sparse linear regression model. Performs the simulation study for both the MCMC and VaNBayes after reading in the above functions.
compile_results.R: Code that produces the graphs and figures.