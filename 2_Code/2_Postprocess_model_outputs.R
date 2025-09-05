library(dplyr)
library(ggplot2)
library(sf)
library(stringi)
library(geojsonsf)
library(ggnewscale)
library(patchwork)
library(here)

txt_size = 10
theme_script <- theme_minimal() +
  theme(axis.text = element_text(colour="black", size = txt_size), 
        axis.ticks = element_line(colour="black"),
        legend.text = element_text(colour="black", size = txt_size),
        legend.title = element_text(colour="black", size = txt_size+1),
        axis.title = element_text(colour="black", size = txt_size+1),
        strip.background = element_rect(colour = "gray", fill = "white"),
        strip.text = element_text(colour = "black", size = txt_size+1))
theme_set(theme_script)

# Load data, source functions
load(here("1_Data", "CHIKV_postprocessing_data.RData"))
source(here("2_Code", "R_functions", "summary_functions.R"))

model_name <- "Main model"

param_est_df <- read.csv(here("3_Output_Data",  "param_eval_estimates.csv"))

if(!file.exists(here("4_Output_Figures"))) dir.create(here("4_Output_Figures"))

# Figures 1D, 2B
plot_incidence_fit(cum_incidence_state_df, cum_incidence_br_df, case_pop_df, state_region_match, age_bin_match)

# Figure 1E
plot_mortality_fit(sim_deaths_df, age_bin_match)

# Figure 2A
plot_sero_fit(sero_data_df)

# Figures 2C, S4
get_plot_lambda(param_est_df, model_name, mean_pop_weights, state_region_match, brazil_map, bahia_constraint=TRUE)

# Figures 2D-E
get_prop_infected(case_pop_df, state_region_match, brazil_map)

# Figure 3A
get_rho_t()

# Figures 3A-D
get_plot_rho_ifr_rr(seed = 346, param_est_df, model_name, case_pop_df, age_bin_match, sim_deaths_df, dengue_rr_df)

