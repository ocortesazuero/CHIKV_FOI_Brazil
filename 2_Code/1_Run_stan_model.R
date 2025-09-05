library(cmdstanr)
library(dplyr)
library(loo)
library(tidybayes)
library(parallel)
library(here)

n_cores <- detectCores() - 2

# Load data
age_bin_match <- read.csv(here("1_Data", "age_bin_match.csv"))
state_region_match <- read.csv(here("1_Data", "state_region_match.csv"))
data <- readRDS(here("1_Data", "data_run_stan_model.rds"))

# Source functions to extract model estimates
source(here("2_Code", "R_functions", "get_model_estimates.R"))

# Compile and run stan model, preferably using a HPC
check_cmdstan_toolchain()
mod <- cmdstan_model(here("2_Code", "Stan", "FOI_model.stan"), pedantic = T)

fit <- mod$sample(data=data, seed=193, 
                  chains=4, parallel_chains=4, 
                  output_dir=here("3_Output_Data"), 
                  iter_sampling=3000, thin=10, refresh=500, iter_warmup=1000)

if(!file.exists(here("3_Output_Data"))) dir.create(here("3_Output_Data"))

# Extract estimates and model convergence evaluation for parameters
var_names <-  c("phi", "log_lambda", "log_ifr", "logit_rho", "logit_rho_t")
cmdstanfit <- as.data.frame(fit$summary(variables = var_names,posterior::default_summary_measures()[1:4], 
                                        quantiles = ~ quantile(., probs = c(0.025,0.975)), 
                                        posterior::default_convergence_measures(), 
                                        .cores=n_cores)) %>%
  rename(med = median,
         lwr = `2.5%`,
         upr = `97.5%`)

nrow(cmdstanfit) == length(which(cmdstanfit$rhat < 1.1))
nrow(cmdstanfit) == length(which(cmdstanfit$ess_bulk > 300))
nrow(cmdstanfit) == length(which(cmdstanfit$ess_tail > 300))

write.csv(cmdstanfit, here("3_Output_Data", "param_eval_estimates.csv"), row.names=F)

# Extract expected values for input data
p_cases <- get_var_sex_loc_time_age(fit, "pLabCases", age_bin_match, state_region_match, n_cores)
p_deaths <- get_var_sex_age(fit, "pDeaths", age_bin_match, n_cores)
p_sero <- get_serofit(fit, n_cores)

# Extract other relevant estimates
tot_infections <- get_var_sex_loc_time_age(fit, "ITot", age_bin_match, state_region_match, n_cores)
immune_pop <- get_var_sex_loc_time_age(fit, "Rg", age_bin_match, state_region_match, n_cores)

reporting_sa <- get_var_sex_age(fit, "mean_rho_sa", age_bin_match, n_cores)
reporting_s <- get_var_sex(fit, "mean_rho_s", n_cores)
reporting_age <- get_var_age(fit, "mean_rho_a", age_bin_match, n_cores)
reporting_time <- get_var_time(fit, "mean_rho_t", n_cores)

reporting_overall <- as.data.frame(fit$summary(variables = c("mean_rho"), "median", quantiles = ~ quantile(., probs = c(0.025,0.975)), .cores=n_cores))
colnames(reporting_overall)[2:4] <- c("med", "lwr", "upr")

ifr_s <- get_var_sex(fit, "mean_ifr_s", n_cores)
ifr_age <- get_var_age(fit, "mean_ifr_a", age_bin_match, n_cores)

ifr_overall <- as.data.frame(fit$summary(variables = c("mean_ifr"), "median", quantiles = ~ quantile(., probs = c(0.025,0.975)), .cores=n_cores))
colnames(ifr_overall)[2:4] <- c("med", "lwr", "upr")

rr_reporting_MvsF_age <- get_var_age(fit, "rr_disease_a", age_bin_match, n_cores)
rr_mortality_MvsF_age <- get_var_age(fit, "rr_death_a", age_bin_match, n_cores)

write.csv(p_cases, here("3_Output_Data", "p_lab_cases.csv"), row.names=F)
write.csv(p_deaths, here("3_Output_Data", "p_deaths.csv"), row.names=F)
write.csv(p_sero, here("3_Output_Data", "seroFit.csv"), row.names=F)
write.csv(tot_infections, here("3_Output_Data", "total_infections.csv"), row.names=F)
write.csv(immune_pop, here("3_Output_Data", "Rg.csv"), row.names=F)
write.csv(reporting_sa, here("3_Output_Data", "mean_rho_sa.csv"), row.names=F)
write.csv(reporting_s, here("3_Output_Data", "mean_rho_s.csv"), row.names=F)
write.csv(reporting_age, here("3_Output_Data", "mean_rho_a.csv"), row.names=F)
write.csv(reporting_time, here("3_Output_Data", "mean_rho_t.csv"), row.names=F)

write.csv(reporting_overall, here("3_Output_Data", "mean_rho.csv"), row.names=F)
write.csv(ifr_s, here("3_Output_Data", "mean_ifr_s.csv"), row.names=F)
write.csv(ifr_age, here("3_Output_Data", "mean_ifr_a.csv"), row.names=F)
write.csv(ifr_overall, here("3_Output_Data", "mean_ifr.csv"), row.names=F)
write.csv(rr_reporting_MvsF_age, here("3_Output_Data", "rr_disease_a.csv"), row.names=F)
write.csv(rr_mortality_MvsF_age, here("3_Output_Data", "rr_death_a.csv"), row.names=F)

# Extract and save chains of several parameters necessary for bootstrapping potential vaccine impact
cmdstanfit <- fit %>%
  gather_draws(logit_rho[Sex,Age_bin_idx]) %>%
  mutate(Sex = ifelse(Sex == 1, "F", "M"),
         Age_bin_idx = Age_bin_idx + 1) %>%
  rename(value = .value,
         variable = .variable)
write.csv(cmdstanfit, here("3_Output_Data", "logit_rho_draws.csv"), row.names=F)

cmdstanfit <- fit %>%
  gather_draws(log_ifr[Sex,Age_bin_idx]) %>%
  mutate(Sex = ifelse(Sex == 1, "F", "M"),
         Age_bin_idx = Age_bin_idx + 1) %>%
  rename(value = .value,
         variable = .variable)
write.csv(cmdstanfit, here("3_Output_Data", "log_ifr_draws.csv"), row.names=F)

cmdstanfit <- fit %>%
  gather_draws(log_lambda[State_code_stan, Year]) %>%
  mutate(Year = Year + 2013) %>%
  rename(value = .value,
         variable = .variable)
write.csv(cmdstanfit, here("3_Output_Data", "log_lambda_draws.csv"), row.names=F)

cmdstanfit <- fit %>%
  gather_draws(phi) %>%
  rename(value = .value,
         variable = .variable)
write.csv(cmdstanfit, here("3_Output_Data", "phi_draws.csv"), row.names=F)

cmdstanfit <- fit %>%
  gather_draws(Rg[Sex, State_code_stan, Year, Age_bin_idx]) %>%
  mutate(Sex = ifelse(Sex == 1, "F", "M"),
         Year = Year + 2013,
         Age_bin_idx = Age_bin_idx + 1) %>%
  rename(value = .value,
         variable = .variable)
write.csv(cmdstanfit,here("3_Output_Data",  "Rg_draws.csv"), row.names=F)

cmdstanfit <- fit %>%
  gather_draws(logit_rho_t[Year]) %>%
  mutate(Year = Year + 2013) %>%
  rename(value = .value,
         variable = .variable)
write.csv(cmdstanfit, here("3_Output_Data", "logit_rho_t_draws.csv"), row.names=F)

ll <- fit$draws('log_lik')
rm(fit)
gc()

# Extract ELPD and DIC
rel_eff <- relative_eff(exp(ll), cores = n_cores)
loo_1 <- loo(ll, cores = n_cores, r_eff = rel_eff)
waicloo <- rbind(as.data.frame(loo_1$estimates), as.data.frame(waic(ll)$estimates))
ll <- posterior::merge_chains(ll)
ll3 <- matrix(ll, nrow=dim(ll)[1], ncol=length(ll)/(dim(ll)[1]))
Dev <- -2*colSums(t(ll3))
DIC <- list(DIC=mean(Dev) + var(Dev)/2, Dbar=mean(Dev), pV=var(Dev)/2)
waicloo[nrow(waicloo)+1,1] <- DIC$DIC
row.names(waicloo)[nrow(waicloo)] <- 'dic'
write.csv(waicloo, here("3_Output_Data", "waicloo.csv"))

















