# Elliot Maceda
# NCSU, May 2024
# Spatial SIR Simulation Study
### Spatial SIR Grid Setup


# only one package necessary, for rdirichlet 
library(gtools)

############################# GRID CREATION ####################################

set.seed(1) # setting this to be constant so we get the same grid each time (later we'll call set.seed(NULL) to undo this)

# lattice dimensions and population
rows = 10
cols = 10
locs = rows*cols

# creating population over a grid
alpha = 1
total_people = 1000
N_vector = as.numeric(round(rdirichlet(1, rep(alpha,locs))*total_people)) # Dirichlet to randomize concentrations of people
N_vector[(N_vector == 0) | (N_vector == 1)] = 2 # make sure there are no holes
N = matrix(N_vector, nrow = rows, byrow = FALSE)

# creating initial infected
I0 = 10
I_vector_index = round(runif(1, 0.5, locs + 0.5))
I_vector = rep(0, locs)
I_vector[I_vector_index] = I0
I = matrix(I_vector, nrow = rows, byrow = FALSE)

# Combining N and I together, calculating S and R.
N = N + I # add the initial infected to the population
S = N - I
R = matrix(0, nrow = rows, ncol = cols)
total_people = sum(N) # recalculate the total number of people (adding I0)


################# CREATING MATRICES OF NEIGHBORING INFECTED ####################


# helps reference each neighbor
col_index = matrix(rep(1:cols,rows), byrow = T, nrow = rows)
row_index = matrix(rep(1:rows,cols), byrow = F, nrow = rows)

# number of infected from each neighbor
left_neighbor_I = matrix(0, nrow = rows, ncol = cols)
right_neighbor_I = matrix(0, nrow = rows, ncol = cols)
top_neighbor_I = matrix(0, nrow = rows, ncol = cols)
bottom_neighbor_I = matrix(0, nrow = rows, ncol = cols)

for (i in 1:rows){
  for (j in 1:cols){
    if (j-1>0){ # left neighbor infected count
      left_neighbor_I[i,j] = I[i,j-1]
    } else{
      left_neighbor_I[i,j] = 0
    }
    
    
    if (j+1<=cols){ # right neighbor infected count
      right_neighbor_I[i,j] = I[i,j+1]
    } else{
      right_neighbor_I[i,j] = 0
    }
    
    
    if (i-1>0){ # top neighbor infected count
      top_neighbor_I[i,j] = I[i-1,j]
    } else{
      top_neighbor_I[i,j] = 0
    }
    
    
    if (i+1<=rows){ # bottom neighbor infected count
      bottom_neighbor_I[i,j] = I[i+1,j]
    } else{
      bottom_neighbor_I[i,j] = 0
    }
  }
}


################# DEFINING VARIABLES FOR THE SIR DIFF EQUATION #################

# (nuisance) parameters for the spatial SIR differential equation 
num_events = 2*total_people # 2 events for each person, "sick" and "recover"
masking_probability = 1
number_reporting_times = 21
delta_t = 0.01

# times we will obtain the data results
reporting_times = round(seq(1, num_events, length.out = number_reporting_times))

# helps with retrieving the data from the correctly chosen times
grid_observation_indices = rep(1:num_events, times = 1, each = rows*cols) %in% reporting_times

# frontloaded variables defined to make the computation a bit faster
location_point_grid = rep(FALSE, rows*cols)  # initializations of containers for interim results
I_grid_count = numeric(rows*cols*num_events) # 
R_grid_count = numeric(rows*cols*num_events) # 

dice_roll = runif(num_events, 0, 1)          # frontloaded realizations of random variables
location_dice_roll = runif(num_events, 0, 1) # 
