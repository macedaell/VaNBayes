adam.R: The optimization algorithm for the neural network.
hetnorm2.R: Includes the initialization, loss function, and predictions for our variational posterior-- the heterogeneous marginal normal model.
SIR_Grid_Setup.R: Helps create the SIR grid
Spatial_SIR_fucntion_large_grids.R: Defines the SIR generating model.
training_distribution_generation.R: Generates data using the training distribution, which is the model's parameter in this case.
zika Brazil data.rdata: Includes data involving the Brazilian States. Originally created from on the Github page: https://github.com/jptrostle/SpatialSIRGPMC
data_analysis.R: The source file. Performs the data analysis by running the above scripts and importing the data generated in "training_distribution_generation".

training_data.csv
observed_data.csv