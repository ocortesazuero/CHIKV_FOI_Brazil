library(dplyr)
library(future.apply)
library(here)
library(parallel)

options(dplyr.summarise.inform = FALSE)

n_workers <- parallel::detectCores() - 2
chunk_size <- n_workers %/% 2
plan(multicore, workers = n_workers)

# We recommend to run this in a HPC

#  Import data, functions get top 5/10 states ---------------------------------------

state_region_match <- read.csv(here("1_Data", "state_region_match.csv"))
age_bin_match <- read.csv(here("1_Data", "age_bin_match.csv"))
life_exp_df <- read.csv(here("1_Data", "life_exp_age_sex.csv"))
sens_characteristics_df <- read.csv(here("1_Data","characteristics_sens_analyses_vacc_impact_forward.csv"))
all_sens_params <- readRDS(here("1_Data","params_sens_analyses_vacc_impact_forward.rda"))

age_group_match <- data.frame(Age_group_idx = 1:5,
                              Age_group = c(">1yo", ">12yo", ">18yo", ">50yo", "18-59"),
                              Age_bin_min = c(2, 4, 5, 12, 5),
                              Age_bin_max = c(16, 16, 16, 16, 13),
                              Age_min = c(1, 12, 18, 50, 18),
                              Age_max = c(90, 90, 90, 90, 59))

source(here("2_Code", "R_functions", "util_functions_vaccine_impact.R"))

top_states <- read.csv(here("1_Data", "case_pop_counts_sex_state_year_age.csv")) %>%
  group_by(State_code_stan) %>%
  summarise(Cases = sum(Cases)) %>%
  arrange(desc(Cases)) %>% slice_head(n=10) %>% select(State_code_stan) %>% as.vector()

# Import population data and aggregate population for calculation of necessary doses
pop_df <- read.csv(here("1_Data", "pop_BR_state_sex_age_2010-2060.csv"), fileEncoding = "latin1") %>%
  mutate(State_code_stan = state_region_match$State_code_stan[match(State, state_region_match$State_ascii)],
         Age = as.numeric(ifelse(Age == "90+", 90, Age))) %>%
  filter(Year %in% c(2024:2029)) 

doses_initial_campaign_df <- bind_rows(lapply(1:5, function(x){
  doses_group_df <- pop_df %>%
    filter(Year == 2024, Age >= age_group_match$Age_min[x], Age <= age_group_match$Age_max[x]) %>%
    group_by(State_code_stan) %>%
    summarise(Pop = sum(Pop)) %>%
    mutate(Age_group_idx = age_group_match$Age_group_idx[x])
  return(doses_group_df)
}))

doses_2dose_campaign_df <- doses_initial_campaign_df %>%
  mutate(Pop <- Pop * 1.9)

doses_catchup_campaign_df <- bind_rows(lapply(1:5, function(x){
  doses_group_df <- pop_df %>%
    filter(Year %in% 2025:2028, Age == age_group_match$Age_min[x]) %>%
    group_by(State_code_stan) %>%
    summarise(Pop_catchup = sum(Pop)) %>%
    mutate(Age_group_idx = age_group_match$Age_group_idx[x])
  return(doses_group_df)
}))

doses_campaign_with_catch_up_df <- doses_initial_campaign_df %>%
  left_join(doses_catchup_campaign_df, by = join_by(State_code_stan, Age_group_idx)) %>%
  mutate(Pop = Pop + Pop_catchup)

pop_df <- pop_df %>%
  mutate(Year = substr(Year, 3, 4)) %>%
  tidyr::pivot_wider(names_from = Year, values_from = Pop, names_prefix = "Pop_") %>%
  mutate(Age_bin_idx = cut(Age, c(age_bin_match$min, 120), right = F, labels = F)) %>%
  select(-c(State_abbr, State))


# Sample relevant parameters from posteriors and prepare objects for bootstrapping --------

set.seed(196)
# FOIs
foi_draws <- read.csv(here("3_Output_Data", "log_lambda_draws.csv")) 

max_draws <- max(foi_draws$.draw)

sample_draws <- sample(max_draws, 1e3, replace = T)
# sample_draws <- sample(max_draws, 5, replace = T)
foi_draws <- foi_draws %>%
  filter(Year >= 2020) %>%
  mutate(lambda = exp(value)) %>%
  group_by(.draw, State_code_stan) %>%
  reframe(foi = sample(lambda, 5, replace = TRUE),
          Year = 25:29) %>% 
  tidyr::pivot_wider(names_from = Year, values_from = foi, names_prefix = "foi_")

foi_draws <- lapply(sample_draws, function(x) foi_draws %>%
                      filter(.draw == x) %>%
                      select(-.draw))

# Population infected at the end of Sept 2024
immun_draws <- read.csv(here("3_Output_Data", "Rg_draws.csv")) %>%
  filter(Year == 2024, .draw %in% sample_draws) %>%
  rename(immunity = value) %>%
  select(-c(Year, .chain, .iteration, variable))

immun_draws <- bind_rows(immun_draws %>%
                           filter(Age_bin_idx == 2) %>%
                           mutate(Age_bin_idx = 1,
                                  immunity = 0),
                         immun_draws) %>%
  arrange(Sex, State_code_stan, Age_bin_idx)

immun_draws <- lapply(sample_draws, function(x) immun_draws %>%
                        filter(.draw == x) %>%
                        select(-.draw))

rho_draws <- read.csv(here("3_Output_Data", "logit_rho_draws.csv")) %>%
  filter(.draw %in% sample_draws) %>%
  mutate(rho = plogis(value)) %>%
  select(-c(.chain, .iteration, variable, value))

rho_t_draws <- read.csv(here("3_Output_Data", "logit_rho_t_draws.csv")) %>%
  filter(.draw %in% sample_draws, Year == 2024) %>%
  mutate(rho_t = plogis(value)) %>%
  select(-c(.chain, .iteration, variable, value, Year))

rho_draws <- rho_draws %>%
  left_join(rho_t_draws, by=".draw") %>%
  mutate(rho= rho*rho_t) %>%
  select(-rho_t)

rho_draws <- bind_rows(rho_draws %>%
                         filter(Age_bin_idx == 2) %>%
                         mutate(Age_bin_idx = 1,
                                rho = 0),
                       rho_draws) %>%
  arrange(Sex, Age_bin_idx, .draw)

rho_draws <- lapply(sample_draws, function(x) rho_draws %>%
                      filter(.draw == x) %>%
                      select(-.draw))

ifr_draws <- read.csv(here("3_Output_Data", "log_ifr_draws.csv")) %>%
  mutate(ifr = exp(value)) %>%
  select(-c(.chain, .iteration, variable, value))

ifr_draws <- bind_rows(ifr_draws %>%
                         filter(Age_bin_idx == 2) %>%
                         mutate(Age_bin_idx = 1,
                                ifr = 0),
                       ifr_draws) %>%
  arrange(Sex, Age_bin_idx, .draw)

ifr_draws <- lapply(sample_draws, function(x) ifr_draws %>%
                      filter(.draw == x) %>%
                      select(-.draw))

gc()

mean_ifr <- read.csv(here("3_Output_Data","mean_ifr_overall.csv"))

ifr_draws_sens_1 <- lapply(ifr_draws, function(ifr_df) ifr_df %>%
                             mutate(ifr = ifelse(Age_bin_idx >=2, ifr/mean_ifr$med*2e-5, ifr)) %>%
                             select(Sex, Age_bin_idx, ifr))

ifr_draws_sens_2 <- lapply(ifr_draws, function(ifr_df) ifr_df %>%
                             mutate(ifr = ifelse(Age_bin_idx >=2, ifr/mean_ifr$med*8e-5, ifr)) %>%
                             select(Sex, Age_bin_idx, ifr))

foi_draws_sens_1 <- lapply(foi_draws, function(foi_df) foi_df %>%
                             mutate(across(starts_with("foi"), ~ .x * 0.5)))
foi_draws_sens_2 <- lapply(foi_draws, function(foi_df) foi_df %>%
                             mutate(across(starts_with("foi"), ~ .x * 1.5)))

foi_draws_sens_3 <- read.csv(here("3_Output_Data","log_lambda_draws.csv")) %>%
  filter(Year >= 2015) %>%
  mutate(lambda = exp(value)) %>%
  group_by(.draw, State_code_stan) %>%
  reframe(foi = sample(lambda, 5, replace = TRUE),
          Year = 25:29) %>% 
  tidyr::pivot_wider(names_from = Year, values_from = foi, names_prefix = "foi_")

foi_draws_sens_3 <- lapply(sample_draws, function(x) foi_draws_sens_3 %>%
                             filter(.draw == x) %>%
                             select(-.draw))

# Base case ------------------------------------------------------

all_states_vacc_effect_df <- sens_summ_effect(immun_draws, foi_draws, rho_draws, ifr_draws, 
                                              append(list(pop_df), all_sens_params[[1]]),
                                              doses_initial_campaign_df,
                                              sens_characteristics_df$type_sens[1], 
                                              sens_characteristics_df$category_sens[1],  
                                              sens_characteristics_df$case_sens[1],  
                                              sens_characteristics_df$bound_sens[1],  
                                              sens_characteristics_df$positive_sens[1],  
                                              sens_characteristics_df$label_sens[1],  
                                              sens_characteristics_df$x_axis_label[1])

for(i in 2:nrow(sens_characteristics_df)){
  if(!(i %in% c(6, 8, 13, 18:22))){
    vacc_effect_df <- sens_summ_effect(immun_draws, foi_draws, rho_draws, ifr_draws, 
                                       append(list(pop_df), all_sens_params[[i]]),
                                       doses_initial_campaign_df,
                                       sens_characteristics_df$type_sens[i], 
                                       sens_characteristics_df$category_sens[i],  
                                       sens_characteristics_df$case_sens[i],  
                                       sens_characteristics_df$bound_sens[i],  
                                       sens_characteristics_df$positive_sens[i],  
                                       sens_characteristics_df$label_sens[i],  
                                       sens_characteristics_df$x_axis_label[i])
  } else if (i == 6){
    vacc_effect_df <- sens_summ_effect(immun_draws, foi_draws, rho_draws, ifr_draws, 
                                       append(list(pop_df), all_sens_params[[i]]),
                                       doses_2dose_campaign_df,
                                       sens_characteristics_df$type_sens[i], 
                                       sens_characteristics_df$category_sens[i],  
                                       sens_characteristics_df$case_sens[i],  
                                       sens_characteristics_df$bound_sens[i],  
                                       sens_characteristics_df$positive_sens[i],  
                                       sens_characteristics_df$label_sens[i],  
                                       sens_characteristics_df$x_axis_label[i])
    
  } else if (i == 8){
    vacc_effect_df <- sens_summ_effect(immun_draws, foi_draws, rho_draws, ifr_draws, 
                                       append(list(pop_df), all_sens_params[[i]]),
                                       doses_campaign_with_catch_up_df,
                                       sens_characteristics_df$type_sens[i], 
                                       sens_characteristics_df$category_sens[i],  
                                       sens_characteristics_df$case_sens[i],  
                                       sens_characteristics_df$bound_sens[i],  
                                       sens_characteristics_df$positive_sens[i],  
                                       sens_characteristics_df$label_sens[i],  
                                       sens_characteristics_df$x_axis_label[i])
    
  } else if (i == 13){
    vacc_effect_df <- sens_summ_effect(immun_draws, foi_draws, rho_draws, ifr_draws, 
                                       append(list(pop_df), all_sens_params[[i]]),
                                       doses_initial_campaign_df,
                                       sens_characteristics_df$type_sens[i], 
                                       sens_characteristics_df$category_sens[i],  
                                       sens_characteristics_df$case_sens[i],  
                                       sens_characteristics_df$bound_sens[i],  
                                       sens_characteristics_df$positive_sens[i],  
                                       sens_characteristics_df$label_sens[i],  
                                       sens_characteristics_df$x_axis_label[i]) %>%
      filter(top_5)
    
  } else if (i == 18){
    vacc_effect_df <- sens_summ_effect(immun_draws, foi_draws, rho_draws, ifr_draws_sens_1, 
                                       append(list(pop_df), all_sens_params[[i]]),
                                       doses_initial_campaign_df,
                                       sens_characteristics_df$type_sens[i], 
                                       sens_characteristics_df$category_sens[i],  
                                       sens_characteristics_df$case_sens[i],  
                                       sens_characteristics_df$bound_sens[i],  
                                       sens_characteristics_df$positive_sens[i],  
                                       sens_characteristics_df$label_sens[i],  
                                       sens_characteristics_df$x_axis_label[i])
  } else if (i == 19){
    vacc_effect_df <- sens_summ_effect(immun_draws, foi_draws, rho_draws, ifr_draws_sens_2, 
                                       append(list(pop_df), all_sens_params[[i]]),
                                       doses_initial_campaign_df,
                                       sens_characteristics_df$type_sens[i], 
                                       sens_characteristics_df$category_sens[i],  
                                       sens_characteristics_df$case_sens[i],  
                                       sens_characteristics_df$bound_sens[i],  
                                       sens_characteristics_df$positive_sens[i],  
                                       sens_characteristics_df$label_sens[i],  
                                       sens_characteristics_df$x_axis_label[i])
  } else if (i == 20){
    vacc_effect_df <- sens_summ_effect(immun_draws, foi_draws_sens_1, rho_draws, ifr_draws, 
                                       append(list(pop_df), all_sens_params[[i]]),
                                       doses_initial_campaign_df,
                                       sens_characteristics_df$type_sens[i], 
                                       sens_characteristics_df$category_sens[i],  
                                       sens_characteristics_df$case_sens[i],  
                                       sens_characteristics_df$bound_sens[i],  
                                       sens_characteristics_df$positive_sens[i],  
                                       sens_characteristics_df$label_sens[i],  
                                       sens_characteristics_df$x_axis_label[i])
  } else if (i == 21){
    vacc_effect_df <- sens_summ_effect(immun_draws, foi_draws_sens_2, rho_draws, ifr_draws, 
                                       append(list(pop_df), all_sens_params[[i]]),
                                       doses_initial_campaign_df,
                                       sens_characteristics_df$type_sens[i], 
                                       sens_characteristics_df$category_sens[i],  
                                       sens_characteristics_df$case_sens[i],  
                                       sens_characteristics_df$bound_sens[i],  
                                       sens_characteristics_df$positive_sens[i],  
                                       sens_characteristics_df$label_sens[i],  
                                       sens_characteristics_df$x_axis_label[i])
  } else if (i == 22){
    vacc_effect_df <- sens_summ_effect(immun_draws, foi_draws_sens_3, rho_draws, ifr_draws, 
                                       append(list(pop_df), all_sens_params[[i]]),
                                       doses_initial_campaign_df,
                                       sens_characteristics_df$type_sens[i], 
                                       sens_characteristics_df$category_sens[i],  
                                       sens_characteristics_df$case_sens[i],  
                                       sens_characteristics_df$bound_sens[i],  
                                       sens_characteristics_df$positive_sens[i],  
                                       sens_characteristics_df$label_sens[i],  
                                       sens_characteristics_df$x_axis_label[i])
  }
  
  all_states_vacc_effect_df <- bind_rows(all_states_vacc_effect_df,
                                         vacc_effect_df)
  print(i)
  
}


write.csv(all_states_vacc_effect_df, here("3_Output_Data", "Vaccine_impact", "all_states_vacc_effect_forward_with_burden_sensitivity_df.csv"), row.names = F)

all_vacc_effect_df <- all_states_vacc_effect_df%>%
  group_by(Type, Category, Case, Bound, Positive, Label, x_axis_lab) %>%
  summarise(across(c(matches("(*)_averted_(lwr|med|upr)$"), Doses), sum)) %>%
  mutate(across(matches("(*)_averted_(lwr|med|upr)$"), ~ .x * 1e4 / Doses, 
                .names = "{paste0(substr(.col, 1, nchar(.col)-3), 'per_10k_doses', substr(.col, nchar(.col)-3, nchar(.col)))}"))

write.csv(all_vacc_effect_df, here("3_Output_Data", "Vaccine_impact", "all_vacc_effect_forward_summary_df.csv"), row.names = F)





