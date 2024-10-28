rm(list=ls())

# packages required
library(geobr) # brazilian state data
library(sf) # for geom_sf
library(dplyr) # for data cleaning
library(ggplot2) # ggplot2


# read all states from 2019 using 'geobr' package
states <- read_state(
  year = 2019, # can change the year if so desired
  showProgress = FALSE
)

head(states)

# rename states so that we don't need to deal with accents when merging the data
states = states %>% 
  mutate(renamed_state = case_when(name_state == "Rondônia" ~ "Rondonia",
                                   name_state == "Amazônas" ~ "Amazonas",
                                   name_state == "Pará" ~ "Para",
                                   name_state == "Amapá" ~ "Amapa",
                                   name_state == "Maranhão" ~ "Maranhao",
                                   name_state == "Piauí" ~ "Piaui",
                                   name_state == "Ceará" ~ "Ceara",
                                   name_state == "Paraíba" ~ "Paraiba",
                                   name_state == "Espírito Santo" ~ "Espirito Santo",
                                   name_state == "São Paulo" ~ "Sao Paulo",
                                   name_state == "Paraná" ~ "Parana",
                                   name_state == "Goiás" ~ "Goias",
                                   .default = as.character(name_state)))

# import our data
load("zika Brazil data.rdata")

# y is 27 x 40... 27 states and 40 time points (weeks)
y <- cbind(zika_data[,, 1], zika_data[,, 2])[, 41:80]
row.names(y) <- states_names

# summarize y somehow (here, I averaged across time)
y_summarized <- apply(y, 1, mean)

# put into a dataframe so we can merge them
y_df <- data.frame(average_infected = unname(y_summarized), 
                   renamed_state = names(y_summarized))

# finally, merge!
states <- merge(states, y_df)

# Remove plot axis for theme_minimal
no_axis <- theme(axis.title=element_blank(),
                 axis.text=element_blank(),
                 axis.ticks=element_blank())


# Plot all Brazilian states without any coloring 
ggplot() +
  geom_sf(data=states, fill="#FFFFFF", color="#000000", linewidth =1, show.legend = FALSE) +
  labs(subtitle="Brazilian States, according to 2019",size=8) +
  theme_void()
  # theme_minimal() + no_axis
# ^ can choose theme_void or theme_minimal by commenting 



# Plot all Brazilian states with coloring according to average_infected
ggplot() +
  geom_sf(data=states, aes(fill=average_infected), color="#000000", linewidth =1) +
  labs(size=8) +
  scale_fill_distiller(palette = "Blues", name="average_infected", limits = c(min(states$average_infected),max(states$average_infected))) +
  theme_void()
# theme_minimal() + no_axis
# ^ can choose theme_void or theme_minimal by commenting


# maybe the palette can be changed? 
