adam.R: The optimization algorithm for the neural network.
hetnorm2.R: Includes the initialization, loss function, and predictions for our variational posterior-- the heterogeneous marginal normal model.
training_distribution_generation.R: Generates data using the training distribution, which is the model's prior in this case. Produces training_data.csv and observed_data.csv.
zika Brazil data.rdata: Includes data involving the Brazilian States. Originally created from on the Github page: https://github.com/jptrostle/SpatialSIRGPMC
data_analysis.R: The source file. Performs the data analysis by running the above scripts and importing the data generated in "training_distribution_generation".
observed_data.csv: The data used in the data analysis, described in Trostle et al. (2024)
map_of_brazilian_states.R: script that was unused in the paper, but could be used to make a map of the brazilian states, colored by average zika virus cases.
