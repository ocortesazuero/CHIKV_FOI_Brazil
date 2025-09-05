transmission_sim <- function(pop_df, immun_df, foi_df, rho_df, ifr_df,
                             prob_disease, prob_mild, prob_severe, prob_chronic,
                             coverage, v_eff_inf, v_eff_dis, 
                             waning = "none", half_life = 5, min_age, max_age, annual_campaign = F){
  
  # #Debug
  # coverage <- 0.4
  # v_eff_inf <- 0.7
  # v_eff_dis <- 0.95
  # waning <- "none"
  # half_life <- 5
  # min_age <- 12
  # annual_campaign <- F
  # immun_df <- immun_draws[[1]]
  # foi_df <- foi_draws[[1]]
  # rho_df <- rho_draws[[1]]
  # ifr_df <- ifr_draws[[1]]
  
  
  # Parameters to compute probability of mild/severe/chronic cases and DALYs
  # prob_disease <- 0.5
  # prob_mild <- 0.88
  # prob_severe <- 1 - prob_mild
  # prob_chronic <- 0.4*prob_severe
  
  weight_mild <- 0.006
  weight_severe <- 0.133
  weight_chronic <- 0.233
  
  dur_mild <- dur_severe <- 6/365
  dur_chronic <- 1
  
  # Useful for waning immunity
  years_since_vacc <- 0
  
  cov_age_df <- data.frame(Age = 0:90) %>%
    mutate(coverage = if_else(Age >= min_age & Age <= max_age, coverage, 0))
  
  sim_df <- pop_df %>%
    left_join(immun_df, by = join_by(State_code_stan, Sex, Age_bin_idx)) %>%
    left_join(foi_df, by = "State_code_stan") %>%
    left_join(rho_df, by = join_by(Sex, Age_bin_idx)) %>%
    left_join(ifr_df, by = join_by(Sex, Age_bin_idx)) %>%
    left_join(life_exp_df, by = join_by(Sex, Age)) %>%
    left_join(cov_age_df, by = "Age") %>%
    mutate(S_nonv_25 = (1-immunity)*(1-coverage),
           S_vacc_25 = (1-immunity)*coverage,
           I_nonv_25 = foi_25*S_nonv_25,
           I_vacc_25 = (1-v_eff_inf)*foi_25*S_vacc_25,
           across(matches("^I_(nonv|vacc)_25$"), ~.x*Pop_25, .names = "{sub('I_', 'Infec_', .col)}"),
           C_obs_nonv_25 = rho*Infec_nonv_25,
           C_obs_vacc_25 = rho*(1-v_eff_dis)/(1-v_eff_inf)*Infec_vacc_25,
           C_all_nonv_25 = prob_disease*Infec_nonv_25,
           C_all_vacc_25 = prob_disease*(1-v_eff_dis)/(1-v_eff_inf)*Infec_vacc_25,
           D_nonv_25 = ifr*Infec_nonv_25,
           D_vacc_25 = ifr*(1-v_eff_dis)/(1-v_eff_inf)*Infec_vacc_25,
           across(matches("^D_(nonv|vacc)_25$"), ~ .x*life_exp, .names = "{sub('D', 'YLL', .col)}")) %>%
    group_by(State_code_stan, Sex) %>%
    arrange(Age, .by_group = TRUE) %>%
    mutate(S_vacc_26 = c(0, S_vacc_25[1:(n()-1)] - I_vacc_25[1:(n()-1)]),
           S_nonv_26 = c(1, S_nonv_25[1:(n()-1)]- I_nonv_25[1:(n()-1)])) %>%
    ungroup()
  
  for(i in 26:29){
    years_since_vacc <- years_since_vacc + 1
    
    switch(waning,
           "none" = NULL,
           "exponential" = {
             v_eff_inf <- v_eff_inf*exp(-years_since_vacc/half_life)
             v_eff_dis <- v_eff_dis*exp(-years_since_vacc/half_life)
           },
           "sigmoid" = {
             v_eff_inf <- v_eff_inf*(1-(years_since_vacc/half_life)**4/(1+(years_since_vacc/half_life)**4))
             v_eff_dis <- v_eff_dis*(1-(years_since_vacc/half_life)**4/(1+(years_since_vacc/half_life)**4))
           })
    
    if(annual_campaign) {
      sim_df <- sim_df %>%
        mutate(across(matches(paste0("^S_vacc_",i,"$")), ~ if_else(Age == min_age, 
                                                                   .x + coverage*get(sub("vacc", "nonv", cur_column())),
                                                                   .x)),
               across(matches(paste0("^S_nonv_",i,"$")), ~ if_else(Age == min_age,
                                                                   (1 - coverage)*.x,
                                                                   .x)))
    }
    
    sim_df <- sim_df %>%
      mutate(across(matches(paste0("^S_nonv_",i,"$")), ~get(paste0("foi_", i))*.x, 
                    .names = "{sub('S', 'I', .col)}"),
             across(matches(paste0("^S_vacc_",i,"$")), ~(1-v_eff_inf)*get(paste0("foi_", i))*.x, 
                    .names = "{sub('S', 'I', .col)}"),
             across(matches(paste0("^I_(nonv|vacc)_",i,"$")), ~.x*get(paste0("Pop_", i)), .names = "{sub('I_', 'Infec_', .col)}"),
             across(matches(paste0("^Infec_nonv_",i,"$")), ~ rho*.x, .names = "{sub('Infec', 'C_obs', .col)}"),
             across(matches(paste0("^Infec_vacc_",i,"$")), ~ (1-v_eff_dis)/(1-v_eff_inf)*rho*.x, 
                    .names = "{sub('Infec', 'C_obs', .col)}"),
             across(matches(paste0("^Infec_nonv_",i,"$")), ~ prob_disease*.x, .names = "{sub('Infec', 'C_all', .col)}"),
             across(matches(paste0("^Infec_vacc_",i,"$")), ~ (1-v_eff_dis)/(1-v_eff_inf)*prob_disease*.x, 
                    .names = "{sub('Infec', 'C_all', .col)}"),
             across(matches(paste0("^Infec_nonv_",i,"$")), ~ ifr*.x, .names = "{sub('Infec', 'D', .col)}"),
             across(matches(paste0("^Infec_vacc_",i,"$")), ~ (1-v_eff_dis)/(1-v_eff_inf)*ifr*.x, 
                    .names = "{sub('Infec', 'D', .col)}"),
             across(matches(paste0("^D_(nonv|vacc)_", i,"$")), ~ .x*life_exp,
                    .names = "{sub('D', 'YLL', .col)}")) %>%
      group_by(State_code_stan, Sex) %>%
      arrange(Age, .by_group = TRUE) %>%
      mutate(across(matches(paste0("^S_nonv_",i,"$")), ~c(1, .x[1:(n()-1)] - get(sub("S_", "I_", cur_column()))[1:(n()-1)]),
                    .names = "{sub(i, i+1, .col)}"),
             across(matches(paste0("^S_vacc_",i,"$")), ~c(0, .x[1:(n()-1)] - get(sub("S_", "I_", cur_column()))[1:(n()-1)]),
                    .names = "{sub(i, i+1, .col)}")) %>%
      ungroup()
    
  }
  
  sim_df <- sim_df %>%
    tidyr::pivot_longer(matches("(Infec|C_all|C_obs|D|YLL)_(vacc|nonv)_(\\d+)"), 
                        names_pattern = "(Infec|C_all|C_obs|D|YLL)_(vacc|nonv)_(\\d+)",
                        names_to = c("Outcome", "Vacc", "Year"),
                        values_to = "Count") %>%
    group_by(State_code_stan, Outcome) %>%
    summarise(Count = sum(Count, na.rm = T)) %>%
    tidyr::pivot_wider(names_from = "Outcome", values_from = "Count")  %>%
    mutate(Mild = prob_mild*C_all,
           Severe = prob_severe*C_all,
           Chronic = prob_chronic*C_all,
           YLD = weight_mild*dur_mild*Mild + weight_severe*dur_severe*Severe + weight_chronic*dur_chronic*Chronic,
           DALYs = YLD + YLL,
           coverage = coverage,
           age_min_vacc = min_age,
           age_max_vacc = max_age)
  
  return(sim_df)
}

transmission_sim_2doses <- function(pop_df, immun_df, foi_df, rho_df, ifr_df,
                                    prob_disease, prob_mild, prob_severe, prob_chronic,
                                    coverage, v_eff_inf, v_eff_dis, 
                                    waning = "none", half_life = 5, min_age, max_age, annual_campaign = F){
  
  # #Debug
  # coverage <- 0.4
  # v_eff_inf <- 0.7
  # v_eff_dis <- 0.95
  # waning <- "none"
  # half_life <- 5
  # min_age <- 12
  # annual_campaign <- F
  # immun_df <- immun_draws[[1]]
  # foi_df <- foi_draws[[1]]
  # rho_df <- rho_draws[[1]]
  # psi_df <- psi_draws[[1]]
  # ifr_df <- ifr_draws[[1]]
  
  
  # Parameters to compute probability of mild/severe/chronic cases and DALYs
  # prob_mild <- 0.88
  # prob_severe <- 1 - prob_mild
  # prob_chronic <- 0.4*prob_severe
  
  weight_mild <- 0.006
  weight_severe <- 0.133
  weight_chronic <- 0.233
  
  dur_mild <- dur_severe <- 6/365
  dur_chronic <- 1
  
  # Useful for waning immunity
  years_since_vacc <- 0
  
  cov_age_df <- data.frame(Age = 0:90) %>%
    mutate(coverage = if_else(Age >= min_age & Age <= max_age, coverage, 0))
  
  sim_df <- pop_df %>%
    left_join(immun_df, by = join_by(State_code_stan, Sex, Age_bin_idx)) %>%
    left_join(foi_df, by = "State_code_stan") %>%
    left_join(rho_df, by = join_by(Sex, Age_bin_idx)) %>%
    left_join(ifr_df, by = join_by(Sex, Age_bin_idx)) %>%
    left_join(life_exp_df, by = join_by(Sex, Age)) %>%
    left_join(cov_age_df, by = "Age") %>%
    mutate(S_nonv_25 = (1-immunity)*(1-coverage),
           S_2vacc_25 = (1-immunity)*coverage*0.9,
           S_1vacc_25 = (1-immunity)*coverage*0.1,
           I_nonv_25 = foi_25*S_nonv_25,
           I_2vacc_25 = (1-v_eff_inf)*foi_25*S_2vacc_25,
           I_1vacc_25 = (1-v_eff_inf*0.5)*foi_25*S_1vacc_25,
           across(matches("^I_(nonv|1vacc|2vacc)_25$"), ~.x*Pop_25, .names = "{sub('I_', 'Infec_', .col)}"),
           C_obs_nonv_25 = rho*Infec_nonv_25,
           C_obs_2vacc_25 = rho*(1-v_eff_dis)/(1-v_eff_inf)*Infec_2vacc_25,
           C_obs_1vacc_25 = rho*(1-v_eff_dis*0.5)/(1-v_eff_inf*0.5)*Infec_1vacc_25,
           C_all_nonv_25 = prob_disease*Infec_nonv_25,
           C_all_2vacc_25 = prob_disease*(1-v_eff_dis)/(1-v_eff_inf)*Infec_2vacc_25,
           C_all_1vacc_25 = prob_disease*(1-v_eff_dis*0.5)/(1-v_eff_inf*0.5)*Infec_1vacc_25,
           D_nonv_25 = ifr*Infec_nonv_25,
           D_2vacc_25 = ifr*(1-v_eff_dis)/(1-v_eff_inf)*Infec_2vacc_25,
           D_1vacc_25 = ifr*(1-v_eff_dis*0.5)/(1-v_eff_inf*0.5)*Infec_2vacc_25,
           across(matches("^D_(nonv|1vacc|2vacc)_25$"), ~ .x*life_exp, .names = "{sub('D', 'YLL', .col)}")) %>%
    group_by(State_code_stan, Sex) %>%
    arrange(Age, .by_group = TRUE) %>%
    mutate(S_2vacc_26 = c(0, S_2vacc_25[1:(n()-1)] - I_2vacc_25[1:(n()-1)]),
           S_1vacc_26 = c(0, S_1vacc_25[1:(n()-1)] - I_1vacc_25[1:(n()-1)]),
           S_nonv_26 = c(1, S_nonv_25[1:(n()-1)]- I_nonv_25[1:(n()-1)])) %>%
    ungroup()
  
  for(i in 26:29){
    years_since_vacc <- years_since_vacc + 1
    
    switch(waning,
           "none" = NULL,
           "exponential" = {
             v_eff_inf <- v_eff_inf*exp(-years_since_vacc/half_life)
             v_eff_dis <- v_eff_dis*exp(-years_since_vacc/half_life)
           },
           "sigmoid" = {
             v_eff_inf <- v_eff_inf*(1-(years_since_vacc/half_life)**4/(1+(years_since_vacc/half_life)**4))
             v_eff_dis <- v_eff_dis*(1-(years_since_vacc/half_life)**4/(1+(years_since_vacc/half_life)**4))
           })
    
    if(annual_campaign) {
      sim_df <- sim_df %>%
        mutate(across(matches(paste0("^S_2vacc_",i,"$")), ~ if_else(Age == min_age, 
                                                                    .x + coverage*0.9*get(sub("2vacc", "nonv", cur_column())),
                                                                    .x)),
               across(matches(paste0("^S_1vacc_",i,"$")), ~ if_else(Age == min_age, 
                                                                    .x + coverage*0.1*get(sub("1vacc", "nonv", cur_column())),
                                                                    .x)),
               across(matches(paste0("^S_nonv_",i,"$")), ~ if_else(Age == min_age,
                                                                   (1 - coverage)*.x,
                                                                   .x)))
    }
    
    sim_df <- sim_df %>%
      mutate(across(matches(paste0("^S_nonv_",i,"$")), ~get(paste0("foi_", i))*.x, 
                    .names = "{sub('S', 'I', .col)}"),
             across(matches(paste0("^S_2vacc_",i,"$")), ~(1-v_eff_inf)*get(paste0("foi_", i))*.x, 
                    .names = "{sub('S', 'I', .col)}"),
             across(matches(paste0("^S_1vacc_",i)), ~(1-v_eff_inf*0.5)*get(paste0("foi_", i))*.x, 
                    .names = "{sub('S', 'I', .col)}"),
             across(matches(paste0("^I_(nonv|1vacc|2vacc)_",i,"$")), ~.x*get(paste0("Pop_", i)), .names = "{sub('I_', 'Infec_', .col)}"),
             across(matches(paste0("^Infec_nonv_",i,"$")), ~ rho*.x, .names = "{sub('Infec', 'C_obs', .col)}"),
             across(matches(paste0("^Infec_2vacc_",i,"$")), ~ (1-v_eff_dis)/(1-v_eff_inf)*rho*.x, 
                    .names = "{sub('Infec', 'C_obs', .col)}"), 
             across(matches(paste0("^Infec_1vacc_",i,"$")), ~ (1-v_eff_dis*0.5)/(1-v_eff_inf*0.5)*rho*.x, 
                    .names = "{sub('Infec', 'C_obs', .col)}"),
             across(matches(paste0("^Infec_nonv_",i,"$")), ~ prob_disease*.x, .names = "{sub('Infec', 'C_all', .col)}"),
             across(matches(paste0("^Infec_2vacc_",i,"$")), ~ (1-v_eff_dis)/(1-v_eff_inf)*prob_disease*.x, 
                    .names = "{sub('Infec', 'C_all', .col)}"), 
             across(matches(paste0("^Infec_1vacc_",i,"$")), ~ (1-v_eff_dis*0.5)/(1-v_eff_inf*0.5)*prob_disease*.x, 
                    .names = "{sub('Infec', 'C_all', .col)}"),
             across(matches(paste0("^Infec_nonv_",i,"$")), ~ ifr*.x, .names = "{sub('Infec', 'D', .col)}"),
             across(matches(paste0("^Infec_2vacc_",i,"$")), ~ (1-v_eff_dis)/(1-v_eff_inf)*ifr*.x, 
                    .names = "{sub('Infec', 'D', .col)}"),
             across(matches(paste0("^Infec_1vacc_",i,"$")), ~ (1-v_eff_dis*0.5)/(1-v_eff_inf*0.5)*ifr*.x, 
                    .names = "{sub('Infec', 'D', .col)}"),
             across(matches(paste0("^D_(nonv|1vacc|2vacc)_", i,"$")), ~ .x*life_exp,
                    .names = "{sub('D', 'YLL', .col)}")) %>%
      group_by(State_code_stan, Sex) %>%
      arrange(Age, .by_group = TRUE) %>%
      mutate(across(matches(paste0("^S_nonv_",i,"$")), ~c(1, .x[1:(n()-1)] - get(sub("S_", "I_", cur_column()))[1:(n()-1)]),
                    .names = "{sub(i, i+1, .col)}"),
             across(matches(paste0("^S_2vacc_",i,"$")), ~c(0, .x[1:(n()-1)] - get(sub("S_", "I_", cur_column()))[1:(n()-1)]),
                    .names = "{sub(i, i+1, .col)}"),
             across(matches(paste0("^S_1vacc_",i,"$")), ~c(0, .x[1:(n()-1)] - get(sub("S_", "I_", cur_column()))[1:(n()-1)]),
                    .names = "{sub(i, i+1, .col)}")) %>%
      ungroup()
    
  }
  
  sim_df <- sim_df %>%
    mutate(across(matches("(Infec|C_all|C_obs|D|YLL)_1vacc_(\\d+)"), ~.x + get(sub("1vacc", "2vacc", cur_column())),
                  .names = "{sub('1vacc', 'vacc', .col)}")) %>%
    tidyr::pivot_longer(matches("(Infec|C_all|C_obs|D|YLL)_(vacc|nonv)_(\\d+)"), 
                        names_pattern = "(Infec|C_all|C_obs|D|YLL)_(vacc|nonv)_(\\d+)",
                        names_to = c("Outcome", "Vacc", "Year"),
                        values_to = "Count") %>%
    group_by(State_code_stan, Outcome) %>%
    summarise(Count = sum(Count, na.rm = T)) %>%
    tidyr::pivot_wider(names_from = "Outcome", values_from = "Count")  %>%
    mutate(Mild = prob_mild*C_all,
           Severe = prob_severe*C_all,
           Chronic = prob_chronic*C_all,
           YLD = weight_mild*dur_mild*Mild + weight_severe*dur_severe*Severe + weight_chronic*dur_chronic*Chronic,
           DALYs = YLD + YLL,
           coverage = coverage,
           age_min_vacc = min_age,
           age_max_vacc = max_age)
  
  return(sim_df)
}

single_scenario_effect <- function(immun_df, foi_df, rho_df, ifr_df,
                                   pop_df, prob_disease, prob_mild, prob_severe, prob_chronic,
                                   coverage, v_eff_inf, v_eff_dis, 
                                   waning, half_life, min_age, max_age, annual_campaign, n_doses) {
  # Wrapper function that runs the transmission_sim function with and with no coverage
  # to evaluate the impact of a vaccination strategy according to sampled parameter estimates,
  # rollout strategies, and vaccine characteristics.
  
  # coverage <- 0.4
  # v_eff_inf <- 0.7
  # v_eff_dis <- 0.95
  # waning <- "none"
  # half_life <- 5
  # min_age <- 12
  # annual_campaign <- F
  # n_doses <- 1
  # immun_df <- immun_draws[[1]]
  # foi_df <- foi_draws[[1]]
  # rho_df <- rho_draws[[1]]
  # ifr_df <- ifr_draws[[1]]
  # prob_disease <- 0.5
  # prob_mild <- 0.88
  # prob_severe <- 0.12
  # prob_chronic <- 0.4*prob_severe
  
  transmission_function <- ifelse(n_doses == 1, transmission_sim, transmission_sim_2doses)
  
  no_vacc_effect_df <- transmission_function(pop_df, immun_df, foi_df, rho_df,ifr_df,
                                             prob_disease, prob_mild, prob_severe, prob_chronic,
                                             coverage = 0, v_eff_inf, v_eff_dis, 
                                             waning, half_life, min_age, max_age, annual_campaign) %>%
    rename_with(.cols = c(C_all, C_obs, Infec, D, Mild, Severe, Chronic, DALYs, YLL, YLD), .fn = ~ paste0("no_vacc_", .x)) %>%
    select(-c(coverage, age_min_vacc, age_max_vacc)) 
  
  vacc_effect_df <- transmission_function(pop_df, immun_df, foi_df, rho_df, ifr_df,
                                          prob_disease, prob_mild, prob_severe, prob_chronic,
                                          coverage, v_eff_inf, v_eff_dis,
                                          waning, half_life, min_age, max_age, annual_campaign) %>%
    left_join(no_vacc_effect_df, by = join_by(State_code_stan)) %>%
    mutate(across(c(C_all, C_obs, Infec, D, Mild, Severe, Chronic, DALYs, YLL, YLD), 
                  ~ get(paste0("no_vacc_", cur_column())) - .x, 
                  .names = "{paste0(.col, '_averted')}"))
  return(vacc_effect_df)
  
}

transform_names <- function(x) {
  replacements <- c(
    "C_averted" = "Cases averted",
    "C_all_averted" = "Cases averted",
    "C_obs_averted" = "Detected cases averted",
    "Infec_averted" = "Infections averted",
    "D_averted" = "Deaths averted",
    "Mild_averted" = "Mild cases averted",
    "Severe_averted" = "Severe cases averted",
    "Chronic_averted" = "Chronic cases averted",
    "DALYs_averted" = "DALYs averted",
    "YLL_averted" = "YLLs averted",
    "YLD_averted" = "YLDs averted",
    "C_averted_per_10k_doses" = "Cases averted per 10,000 doses",
    "C_all_averted_per_10k_doses" = "All cases averted per 10,000 doses",
    "C_obs_averted_per_10k_doses" = "Detected cases averted per 10,000 doses",
    "Infec_averted_per_10k_doses" = "Infections averted per 10,000 doses",
    "D_averted_per_10k_doses" = "Deaths averted per 10,000 doses",
    "Mild_averted_per_10k_doses" = "Mild cases averted per 10,000 doses",
    "Severe_averted_per_10k_doses" = "Severe cases averted per 10,000 doses",
    "Chronic_averted_per_10k_doses" = "Chronic cases averted per 10,000 doses",
    "DALYs_averted_per_10k_doses" = "DALYs averted per 10,000 doses",
    "YLL_averted_per_10k_doses" = "YLLs averted per 10,000 doses",
    "YLD_averted_per_10k_doses" = "YLDs averted per 10,000 doses"
  )
  result <- replacements[x]
  return(result)
}

sens_summ_effect <- function(immun_draws, foi_draws, rho_draws, ifr_draws, param_list, doses_campaign_df, type_sens, 
                             category_sens, case_sens, bound_sens, positive_sens, label_sens, x_axis_label){
  sens_effect_df <- bind_rows(future_mapply(single_scenario_effect,
                                            immun_draws,
                                            foi_draws,
                                            rho_draws,
                                            ifr_draws,
                                            MoreArgs = param_list,
                                            SIMPLIFY = FALSE, future.seed = F, future.chunk.size = chunk_size)) #%>%
  unwanted_cols <- c("State_code_stan", "coverage", "age_min_vacc", "age_max_vacc", "Age_bin_min", "Age_bin_max", "Age_group")
  cols_to_sum <- setdiff(names(sens_effect_df), unwanted_cols)
  sens_effect_df <- sens_effect_df %>%
    left_join(age_group_match, by = c("age_min_vacc"="Age_min", "age_max_vacc"="Age_max")) %>%
    group_by(State_code_stan, Age_group_idx) %>%
    summarise(across(all_of(cols_to_sum), ~ quantile(.x, probs = c(0.5)), .names = "{paste0(.col, '_med')}"),
              across(all_of(cols_to_sum), ~ quantile(.x, probs = c(0.025)), .names = "{paste0(.col, '_lwr')}"),
              across(all_of(cols_to_sum), ~ quantile(.x, probs = c(0.975)), .names = "{paste0(.col, '_upr')}")) %>%
    left_join(doses_campaign_df, by = join_by(State_code_stan, Age_group_idx)) %>%
    mutate(Doses = Pop * param_list$coverage,
           Type = type_sens, Category = category_sens, Case = case_sens, Bound = bound_sens, Positive = positive_sens, 
           Label = label_sens, x_axis_lab = x_axis_label,
           top_5 = State_code_stan %in% top_states[[1]][1:5]) %>%
    ungroup()
  return(sens_effect_df)
}
