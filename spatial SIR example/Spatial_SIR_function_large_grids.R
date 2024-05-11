
# Elliot Maceda
# NCSU, May 2024
# Spatial SIR Simulation Study
### Spatial SIR Stochastic Differential Equation Function

############################ PACKAGES REQUIRED #################################

# no packages are actually necessary - unless you want to make plots:
library(ggplot2)     # for the aggregate time series plot
library(tidyverse)   # 

library(countcolors) # for the RGB lattice plots


########################## SIMULATION FUNCTION #################################
  

# inputs: 
# SIR Stochastic Diff Equation parameters: beta -- local infection rate, phi -- spatial infection rate, eta -- recovery rate
# SIR Grid lattice: N, S, I, R, rows, columns, row_index, col_index, total_people, matrices of infected neighbors
# SIR Stochastic Diff Equation Nuisance parameters: num_events, masking_prob, number_reporting_times, grid_observation_indices
# Frontloaded variables: dice_roll, location_dice_roll, location_point_grid, I_grid_count, R_grid_count
# options to graph or not

Spatial_SIR_simulation = function(beta, phi, eta,
                                  N,S,I,R,rows,cols,row_index, col_index, total_people,
                                  left_neighbor_I,right_neighbor_I,top_neighbor_I,bottom_neighbor_I,
                                  num_events,masking_probability,number_reporting_times,reporting_times,grid_observation_indices,
                                  dice_roll,location_dice_roll, delta_t,
                                  location_point_grid,I_grid_count,R_grid_count,
                                  aggregate_ts_plot,lattice_snapshots){
  
  
  # initialize aggregate time series plot
  if (aggregate_ts_plot){
    S_total_count = numeric(num_events)
    I_total_count = numeric(num_events)
    R_total_count = numeric(num_events)
  }
  
  # initialize RGB lattice plots
  if (lattice_snapshots){
    lattice_snapshot_divisor = array(c(N,N,N), dim = c(rows,cols,3))
    lattice_snapshot = vector(mode = "list", length = num_events)
  }
  
  # helps with main computation
  beta_over_N = beta/N
  phi_over_N = phi/N
  
  # run through the events of individuals becoming infected or recovering
  event = 0
  while (event  < num_events){ 
    event = event + 1
    if (sum(I) != 0){ # if there are still infected people...
      
      # calculate probability of events across the entire grid (infection or recovery)
      prob_newly_infected = S*I*beta_over_N*delta_t + delta_t*S*(left_neighbor_I + right_neighbor_I + top_neighbor_I + bottom_neighbor_I)*phi_over_N
      prob_newly_recovered = eta*I*delta_t
      
      # calculate how much time until this next event using geometric distribution with prob_nothing_happens
      prob_nothing_happens = 1 - sum(prob_newly_infected) - sum(prob_newly_recovered)
      if (prob_nothing_happens < 0){prob_nothing_happens=0} # ensure "prob_nothing_happens" is a valid probability
      time_until_next_event = rgeom(1,1-prob_nothing_happens)
      if ((time_until_next_event > 0)){ # if we need to wait, save the current results across the waiting time
        
        I_grid_count[(rows*cols*(event-1)+1):(rows*cols*(event + time_until_next_event - 1))] = as.vector(I)
        R_grid_count[(rows*cols*(event-1)+1):(rows*cols*(event + time_until_next_event - 1))] = as.vector(R)
        
        if (aggregate_ts_plot){
          S_total_count[event:(event + time_until_next_event - 1)] = sum(S)
          I_total_count[event:(event + time_until_next_event - 1)] = sum(I)
          R_total_count[event:(event + time_until_next_event - 1)] = sum(R)
        }
        
        
        if (lattice_snapshots){
          for (i in event:(event + time_until_next_event - 1)){
            lattice_snapshot[[i]] = (array(c(I,R,S), dim = c(rows,cols,3))/lattice_snapshot_divisor)
          }
        }
        
        # fast forward the clock until the event occurs
        event = event + time_until_next_event
        num_events = num_events + time_until_next_event
        
        # need to generate more realizations for the time lost (can probably implement this in a better way)
        dice_roll = append(dice_roll, runif(time_until_next_event))
        location_dice_roll = append(location_dice_roll, runif(time_until_next_event))
      }
      
      # determine which event will happen now-- a new infection or new recovery
      total_prob_value = sum(prob_newly_infected) + sum(prob_newly_recovered)
      prob_newly_infected_STANDARDIZED = prob_newly_infected/total_prob_value
      prob_newly_recovered_STANDARDIZED = prob_newly_recovered/total_prob_value

      if (dice_roll[event] < sum(prob_newly_infected_STANDARDIZED)){ # in the event of a newly infected...
        
        # determine where this new infection appears
        location_prob = cumsum(prob_newly_infected_STANDARDIZED)/max(cumsum(prob_newly_infected_STANDARDIZED))
        location = which.max(location_dice_roll[event] < location_prob)
        location_point = location_point_grid
        location_point[location]= TRUE
        location_roll = matrix(location_point, nrow = rows, ncol = cols)
        
        # update the lattice grid
        S = S - location_roll
        I = I + location_roll
        
      } else{ # in the event of a recovered individual...
        
        # determine where this new recovery appears
        location_prob = cumsum(prob_newly_recovered_STANDARDIZED)/max(cumsum(prob_newly_recovered_STANDARDIZED))
        location = which.max(location_dice_roll[event] < location_prob)
        location_point = location_point_grid
        location_point[location]= TRUE
        location_roll = matrix(location_point, nrow = rows, ncol = cols)
        
        # update the lattice grid
        I = I - location_roll
        R = R + location_roll
      }
      
      # update the matrices of neighbor's infected counts
      if (col_index[location_roll]+1 <= cols){
        left_neighbor_I[row_index[location_roll],col_index[location_roll]+1] = I[row_index[location_roll],col_index[location_roll]]
      }
      if (0 < col_index[location_roll]-1){
        right_neighbor_I[row_index[location_roll],col_index[location_roll]-1] = I[row_index[location_roll],col_index[location_roll]]
      }
      if (row_index[location_roll]+1 <= rows){
        top_neighbor_I[row_index[location_roll]+1,col_index[location_roll]] = I[row_index[location_roll],col_index[location_roll]]
      }
      if (0 < row_index[location_roll]-1){
        bottom_neighbor_I[row_index[location_roll]-1,col_index[location_roll]] = I[row_index[location_roll],col_index[location_roll]]
      }
      
      # save the data for the current time point
      I_grid_count[(rows*cols*(event-1)+1):(rows*cols*(event))] = as.vector(I)
      R_grid_count[(rows*cols*(event-1)+1):(rows*cols*(event))] = as.vector(R)
      
      # save the data for the aggregate time series plot
      if (aggregate_ts_plot){
        S_total_count[event] = sum(S)
        I_total_count[event] = sum(I)
        R_total_count[event] = sum(R)
      }
      
      # save the data for the RGB lattice plot
      if (lattice_snapshots){
        lattice_snapshot[[event]] = (array(c(I,R,S), dim = c(rows,cols,3))/lattice_snapshot_divisor)
      }
      
      
    }else{ # if the number of infected hits 0 early, save the results over the entire time window and stop the simulation
      
      I_grid_count[(rows*cols*(event-1)+1):(rows*cols*(num_events))] = rep(as.vector(I), num_events - event + 1)
      R_grid_count[(rows*cols*(event-1)+1):(rows*cols*(num_events))] = rep(as.vector(R), num_events - event + 1)
      if (aggregate_ts_plot){
        S_total_count[event:num_events] = S_total_count[event-1]
        I_total_count[event:num_events] = I_total_count[event-1]
        R_total_count[event:num_events] = R_total_count[event-1]
      }
      
      break
    }
  }
  
  
  
  # Simulation output! If masking probability < 1, we may only see a proportion of the infection counts
  Observed_I_grid_count = rbinom(rows*cols*number_reporting_times, I_grid_count[grid_observation_indices], masking_probability)
  Observed_R_grid_count = rbinom(rows*cols*number_reporting_times, R_grid_count[grid_observation_indices], masking_probability)
  
  ########################## DATAVIZ OUTPUT ######################################
  
  
  ### Aggregate Dataviz Output
  if (aggregate_ts_plot){
    
    # organizing the plot's dataframe
    summarizing_df = data.frame(Observed_I_grid_count, Observed_R_grid_count, grouping_vector = as.character(rep(1:number_reporting_times, times = 1, each = rows*cols)), df_reporting_times = rep(reporting_times, times = 1, each = rows*cols))
    summarizing_df2 = summarizing_df%>% 
      group_by(grouping_vector) %>% 
      summarize(Observed_I_total = sum(Observed_I_grid_count),
                Observed_R_total = sum(Observed_R_grid_count)) %>% 
      ungroup() %>% 
      mutate(reporting_times = as.numeric(grouping_vector)) %>% 
      arrange(reporting_times)
    
    # print the plot
    print(ggplot() +
      geom_line(aes(x = (1:num_events), y = (total_people - I_total_count - R_total_count)), color = "blue") +
      geom_line(aes(x = (1:num_events), y = I_total_count), color = "red") +
      geom_line(aes(x = (1:num_events), y = R_total_count), color = "green")  + 
      geom_point(aes(x = reporting_times, y = (total_people - summarizing_df2$Observed_I_total - summarizing_df2$Observed_R_total)), color = "blue") +
      geom_point(aes(x = reporting_times, y = summarizing_df2$Observed_I_total), color = "red") +
      geom_point(aes(x = reporting_times, y = summarizing_df2$Observed_R_total), color = "green") +
      xlab("Time") + ylab("Total Count") + theme_minimal())
  }
  
  
  
  
  # RGB Lattice Plot Output
  if (lattice_snapshots){
    
    last_event = event - 1
    
    # start of the SIR grid
    op <- graphics::par(mar = c(0, 0, 2, 0))
    asp <- dim(lattice_snapshot[[1]])[1]/dim(lattice_snapshot[[1]])[2]
    graphics::plot(0:1, 0:1, type = "n", ann = F, axes = F, asp = asp)
    graphics::rasterImage(lattice_snapshot[[1]], 0, 0, 1, 1, interpolate = FALSE)
    title("Time = 1")
    
    # after 1% of the total time has passed
    op <- graphics::par(mar = c(0, 0, 2, 0))
    asp <- dim(lattice_snapshot[[round(quantile(1:last_event, .01))]])[1]/dim(lattice_snapshot[[round(quantile(1:last_event, .01))]])[2]
    graphics::plot(0:1, 0:1, type = "n", ann = F, axes = F, asp = asp)
    graphics::rasterImage(lattice_snapshot[[round(quantile(1:last_event, .01))]], 0, 0, 1, 1, interpolate = FALSE)
    title(paste0("Time = ",round(quantile(1:last_event, .01))))
    
    # after 5% of the total time has passed
    op <- graphics::par(mar = c(0, 0, 2, 0))
    asp <- dim(lattice_snapshot[[round(quantile(1:last_event, .05))]])[1]/dim(lattice_snapshot[[round(quantile(1:last_event, .05))]])[2]
    graphics::plot(0:1, 0:1, type = "n", ann = F, axes = F, asp = asp)
    graphics::rasterImage(lattice_snapshot[[round(quantile(1:last_event, .05))]], 0, 0, 1, 1, interpolate = FALSE)
    title(paste0("Time = ",round(quantile(1:last_event, .05))))
    
    # after 15% of the total time has passed
    op <- graphics::par(mar = c(0, 0, 2, 0))
    asp <- dim(lattice_snapshot[[round(quantile(1:last_event, .15))]])[1]/dim(lattice_snapshot[[round(quantile(1:last_event, .15))]])[2]
    graphics::plot(0:1, 0:1, type = "n", ann = F, axes = F, asp = asp)
    graphics::rasterImage(lattice_snapshot[[round(quantile(1:last_event, .15))]], 0, 0, 1, 1, interpolate = FALSE)
    title(paste0("Time = ",round(quantile(1:last_event, .15))))
    
    # after 25% of the total time has passed
    op <- graphics::par(mar = c(0, 0, 2, 0))
    asp <- dim(lattice_snapshot[[round(quantile(1:last_event, .25))]])[1]/dim(lattice_snapshot[[round(quantile(1:last_event, .25))]])[2]
    graphics::plot(0:1, 0:1, type = "n", ann = F, axes = F, asp = asp)
    graphics::rasterImage(lattice_snapshot[[round(quantile(1:last_event, .25))]], 0, 0, 1, 1, interpolate = FALSE)
    title(paste0("Time = ",round(quantile(1:last_event, .25))))
    
    # after 50% of the total time has passed
    op <- graphics::par(mar = c(0, 0, 2, 0))
    asp <- dim(lattice_snapshot[[round(quantile(1:last_event, .5))]])[1]/dim(lattice_snapshot[[round(quantile(1:last_event, .5))]])[2]
    graphics::plot(0:1, 0:1, type = "n", ann = F, axes = F, asp = asp)
    graphics::rasterImage(lattice_snapshot[[round(quantile(1:last_event, .5))]], 0, 0, 1, 1, interpolate = FALSE)
    title(paste0("Time = ",round(quantile(1:last_event, .5))))
    
    # after 75% of the total time has passed
    op <- graphics::par(mar = c(0, 0, 2, 0))
    asp <- dim(lattice_snapshot[[round(quantile(1:last_event, .75))]])[1]/dim(lattice_snapshot[[round(quantile(1:last_event, .75))]])[2]
    graphics::plot(0:1, 0:1, type = "n", ann = F, axes = F, asp = asp)
    graphics::rasterImage(lattice_snapshot[[round(quantile(1:last_event, .75))]], 0, 0, 1, 1, interpolate = FALSE)
    title(paste0("Time = ",round(quantile(1:last_event, .75))))
    
    # SIR grid at the end of the time period
    op <- graphics::par(mar = c(0, 0, 2, 0))
    asp <- dim(lattice_snapshot[[last_event]])[1]/dim(lattice_snapshot[[last_event]])[2]
    graphics::plot(0:1, 0:1, type = "n", ann = F, axes = F, asp = asp)
    graphics::rasterImage(lattice_snapshot[[last_event]], 0, 0, 1, 1, interpolate = FALSE)
    title(paste0("Time = ",last_event))
    
    # What the total population layout looks like
    op <- graphics::par(mar = c(0, 0, 2, 0))
    asp <- dim(lattice_snapshot_divisor/max(lattice_snapshot_divisor))[1]/dim(lattice_snapshot_divisor/max(lattice_snapshot_divisor))[2]
    graphics::plot(0:1, 0:1, type = "n", ann = F, axes = F, asp = asp)
    graphics::rasterImage(lattice_snapshot_divisor/max(lattice_snapshot_divisor), 0, 0, 1, 1, interpolate = FALSE)
    title(paste0("Population"))
  }
  
  # end of function, return counts of infected and recovered
  return(response = c(Observed_I_grid_count, Observed_R_grid_count))
  
}


# # for testing:
# 
# true_eta = .1
# true_beta = .1
# true_phi = .1
# delta_t = 0.01
# 
# true_param_test = Spatial_SIR_simulation(true_beta,true_phi,true_eta,
#                                          N,S,I,R,rows,cols,row_index, col_index, total_people,
#                                          left_neighbor_I,right_neighbor_I,top_neighbor_I,bottom_neighbor_I,
#                                          num_events,masking_probability,number_reporting_times,reporting_times,grid_observation_indices,
#                                          dice_roll,location_dice_roll,delta_t,
#                                          location_point_grid,I_grid_count,R_grid_count,
#                                          aggregate_ts_plot = T,lattice_snapshots = T)
# 
# 
# true_param_test
# 
