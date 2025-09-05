library(dplyr)
library(future.apply)
library(here)
library(parallel)

options(dplyr.summarise.inform = FALSE)

# We recommend to run this in a HPC
# n_workers <- chunk_size <- 1
n_workers <- parallel::detectCores() - 2
chunk_size <- n_workers %/% 2
plan(multicore, workers = n_workers)

#  Import data, functions get top 5/10 states ---------------------------------------

state_region_match <- read.csv(here("1_Data", "state_region_match.csv"))
age_bin_match <- read.csv(here("1_Data", "age_bin_match.csv"))
life_exp_df <- read.csv(here("1_Data", "life_exp_age_sex.csv"))

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

# Transmission model and vaccine impact -----------------------------------

# Conditions for Base Case
cov <- 0.4
ve_i <- 0.7
ve_d <- 0.95

bc_params <- list(pop_df = pop_df,
                  prob_disease = 0.5,
                  prob_mild = 0.88,
                  prob_severe = 0.12,
                  prob_chronic = 0.4*0.12,
                  coverage = cov, 
                  v_eff_inf = ve_i,
                  v_eff_dis = ve_d, 
                  waning = "none", 
                  half_life = 3,
                  min_age = 12,
                  max_age = 90,
                  annual_campaign = F,
                  n_doses = 1)

sens_min_age_50_params <- bc_params
sens_min_age_50_params$min_age <- 50

# Bootstrap for base case scenario
all_states_vacc_effect_bc_df <- bind_rows(future_mapply(single_scenario_effect,
                                                        immun_draws,
                                                        foi_draws,
                                                        rho_draws,
                                                        ifr_draws,
                                                        MoreArgs = bc_params,
                                                        SIMPLIFY = FALSE, future.seed = NULL, future.chunk.size = chunk_size))

# Bootstrap for alternative scenario Fig 4
all_states_vacc_effect_50_df <- bind_rows(future_mapply(single_scenario_effect,
                                                        immun_draws,
                                                        foi_draws,
                                                        rho_draws,
                                                        ifr_draws,
                                                        MoreArgs = sens_min_age_50_params,
                                                        SIMPLIFY = FALSE, future.seed = NULL, future.chunk.size = chunk_size))

cols_burden <- colnames(all_states_vacc_effect_bc_df)[!(colnames(all_states_vacc_effect_bc_df) %in% c("State_code_stan",
                                                                                                      "age_min_vacc",
                                                                                                      "age_max_vacc"))]

# Summarize results per state
all_states_vacc_effect_df <- bind_rows(all_states_vacc_effect_bc_df,
                                       all_states_vacc_effect_50_df) %>%
  left_join(age_group_match, by = c("age_min_vacc"="Age_min", "age_max_vacc"="Age_max")) %>%
  group_by(State_code_stan, Age_group_idx) %>%
  summarise(across(all_of(cols_burden), ~ quantile(.x, probs = c(0.5)), .names = "{paste0(.col, '_med')}"),
            across(all_of(cols_burden), ~ quantile(.x, probs = c(0.025)), .names = "{paste0(.col, '_lwr')}"),
            across(all_of(cols_burden), ~ quantile(.x, probs = c(0.975)), .names = "{paste0(.col, '_upr')}")) %>%
  left_join(doses_initial_campaign_df, by = join_by(State_code_stan, Age_group_idx)) %>%
  mutate(Doses = Pop * cov,
         top_5 = State_code_stan %in% top_states[[1]][1:5])

# Summarize results nationally
brazil_impact_df <- all_states_vacc_effect_df %>%
  group_by(Age_group_idx) %>%
  summarise(across(c(matches("(*)_averted_(lwr|med|upr)$"), Doses), sum)) %>%
  mutate(across(matches("(*)_averted_(lwr|med|upr)$"), ~ .x * 1e4 / Doses, 
                .names = "{paste0(substr(.col, 1, nchar(.col)-3), 'per_10k_doses', substr(.col, nchar(.col)-3, nchar(.col)))}")) %>%
  mutate(Subset = "Brazil")

# Top 5 states in terms of cases
top_5_states_impact_df <- all_states_vacc_effect_df %>%
  filter(top_5) %>%
  group_by(Age_group_idx) %>%
  summarise(across(c(matches("(*)_averted_(lwr|med|upr)$"), Doses), sum)) %>%
  mutate(across(matches("(*)_averted_(lwr|med|upr)$"), ~ .x * 1e4 / Doses, 
                .names = "{paste0(substr(.col, 1, nchar(.col)-3), 'per_10k_doses', substr(.col, nchar(.col)-3, nchar(.col)))}")) %>%
  mutate(Subset = "Top 5 States")

# Tidy data frames and get them ready for plotting
vacc_impact_main_df <- bind_rows(brazil_impact_df, top_5_states_impact_df) %>%
  mutate(Age_group= age_group_match$Age_group[match(Age_group_idx, age_group_match$Age_group_idx)]) %>%
  mutate(x = paste(Age_group, Subset, sep = " - ")) %>%
  arrange(Age_group)

vacc_impact_main_plot <- bind_rows(vacc_impact_main_df %>%
                                     tidyr::pivot_longer(cols = matches("(lwr|med|upr)$"), 
                                                         names_to = c("name", ".value"),
                                                         names_pattern = "(.+)_(med|lwr|upr)") %>%
                                     mutate(name = transform_names(name)) %>%
                                     select(-Doses),
                                   vacc_impact_main_df %>%
                                     rename(med = Doses) %>%
                                     mutate(lwr = NA, upr = NA, name = "Doses required") %>%
                                     select(-matches("_(lwr|med|upr)$")))

dir.create(here("3_Output_Data", "Vaccine_impact"))
write.csv(all_states_vacc_effect_df, here("3_Output_Data", "Vaccine_impact", "all_states_vacc_effect_forward_BC_df.csv"), row.names = F)
write.csv(vacc_impact_main_df, here("3_Output_Data", "Vaccine_impact", "vacc_impact_main_observed_cases_df.csv"), row.names = F)
write.csv(vacc_impact_main_plot, here("3_Output_Data", "Vaccine_impact", "vacc_impact_main_plot_observed_cases.csv"), row.names = F)
