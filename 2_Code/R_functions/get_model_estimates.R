get_var_sex_loc_time_age <- function(fit, var, age_bin_match, state_region_match, n_cores){
  cmdstanfit <- as.data.frame(fit$summary(variables = var, "median", quantiles = ~ quantile(., probs = c(0.025,0.975)), .cores=n_cores))
  colnames(cmdstanfit)[2:4] <- c("med", "lwr", "upr")
  
  cmdstanfit[,c("Sex", "State_code_stan", "Year", "Age_bin_idx")] <- sapply(1:4, function(x) readr::parse_number(stringr::str_split_i(cmdstanfit$variable, ",", x)))
  
  cmdstanfit <- cmdstanfit %>%
    mutate(Sex = ifelse(Sex == 1, "F", "M"),
           Age_bin_idx = Age_bin_idx + 1,
           Year = Year + 2013,
           State = state_region_match$state[match(State_code_stan, state_region_match$State_code_stan)],
           Age = age_bin_match$Age[match(Age_bin_idx, age_bin_match$Age_bin_idx)])
  return(cmdstanfit)
}

get_var_sex_age <- function(fit, var, age_bin_match, n_cores){
  cmdstanfit <- as.data.frame(fit$summary(variables = var, "median", quantiles = ~ quantile(., probs = c(0.025,0.975)), .cores=n_cores)) 
  colnames(cmdstanfit)[2:4] <- c("med", "lwr", "upr")
  
  cmdstanfit[,c("Sex", "Age_bin_idx")] <- sapply(1:2, function(x) readr::parse_number(stringr::str_split_i(cmdstanfit$variable, ",", x)))
  
  cmdstanfit <- cmdstanfit %>%
    mutate(Sex = ifelse(Sex == 1, "F", "M"),
           Age_bin_idx = Age_bin_idx + 1,
           Age = age_bin_match$Age[match(Age_bin_idx, age_bin_match$Age_bin_idx)])
  return(cmdstanfit)
}

get_serofit <- function(fit, n_cores){
  cmdstanfit <- as.data.frame(fit$summary(variables = c("meanSero"), "median", quantiles = ~ quantile(., probs = c(0.025,0.975)), .cores=n_cores)) %>%
    rename(med = median,
           lwr = `2.5%`,
           upr = `97.5%`)
  return(cmdstanfit)
}

get_var_age <- function(fit, var, age_bin_match, n_cores){
  cmdstanfit <- as.data.frame(fit$summary(variables = var, "median", quantiles = ~ quantile(., probs = c(0.025,0.975)), .cores=n_cores))
  colnames(cmdstanfit)[2:4] <- c("med", "lwr", "upr")
  
  cmdstanfit[,c("Age_bin_idx")] <- readr::parse_number(stringr::str_split_i(cmdstanfit$variable, ",", 1))
  
  cmdstanfit <- cmdstanfit %>%
    mutate(Age_bin_idx = Age_bin_idx + 1,
           Age = age_bin_match$Age[match(Age_bin_idx, age_bin_match$Age_bin_idx)]) %>%
    select(c(Age_bin_idx, Age, med:upr))
  return(cmdstanfit)
}

get_var_sex <- function(fit, var, n_cores){
  cmdstanfit <- as.data.frame(fit$summary(variables = c("mean_ifr_s"), "median", quantiles = ~ quantile(., probs = c(0.025,0.975)), .cores=n_cores))
  colnames(cmdstanfit)[2:4] <- c("med", "lwr", "upr")
  
  cmdstanfit[,c("Sex")] <- readr::parse_number(stringr::str_split_i(cmdstanfit$variable, ",", 1))
  
  cmdstanfit <- cmdstanfit %>%
    mutate(Sex = ifelse(Sex == 1, "F", "M")) %>%
    select(c(Sex, med:upr))
  return(cmdstanfit)
}

get_var_time <- function(fit, var, n_cores){
  cmdstanfit <- as.data.frame(fit$summary(variables = c("mean_rho_t"), "median", quantiles = ~ quantile(., probs = c(0.025,0.975)), .cores=n_cores))
  colnames(cmdstanfit)[2:4] <- c("med", "lwr", "upr")
  
  cmdstanfit[,c("Year")] <- readr::parse_number(stringr::str_split_i(cmdstanfit$variable, ",", 1))
  
  cmdstanfit <- cmdstanfit %>%
    mutate(Year = Year + 2013) %>%
    select(c(Year, med:upr))
  return(cmdstanfit)
  
}