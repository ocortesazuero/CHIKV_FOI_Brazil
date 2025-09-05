plot_incidence_fit <- function(cum_incidence_state_df, cum_incidence_br_df, case_pop_df, state_region_match, age_bin_match){
  p_cases <- read.csv(here("3_Output_Data", "p_lab_cases.csv")) %>% 
    mutate(Type = "Model") %>%
    select(-variable)
  
  fill_vals <- c("black")
  
  if(min(p_cases$Age_bin_idx) == 2){
    p_cases <- bind_rows(p_cases,
                         case_pop_df %>%
                           filter(Age_bin_idx == 1) %>%
                           rename(med = Cases) %>%
                           select(-c(Sex_stan, Year_stan, State, Pop)) %>%
                           mutate(Type = "Data")) %>%
      mutate(Type = factor(Type, levels = c("Model", "Data")))
    
    fill_vals <- c(fill_vals, "white")
  }
  
  p_cases <- p_cases %>%
    left_join((case_pop_df %>% select(c(Sex, State_code_stan, Year, Age_bin_idx, Cases, Pop))), by = join_by(Sex, State_code_stan, Year, Age_bin_idx)) %>%
    mutate(Incidence = Cases*1e5/11/Pop,
           across(c(med, lwr, upr), ~.x*1e5/11/Pop, .names = "Incidence_{.col}"))
  
  colorpal <- c("#00b9a7", "#ffce00")
  
  dir.create(here("4_Output_Figures", "Incidence"))
  
  for(i in 1:27){
    state_df <- p_cases %>%
      filter(State_code_stan == i)
    
    p <- ggplot(state_df, aes(x=Age_bin_idx)) +
      facet_wrap(. ~ Year) +
      geom_bar(aes(y=Incidence, fill=Sex), stat="identity", position="dodge") +
      scale_fill_manual(values=colorpal) +
      guides(fill = guide_legend(nrow = 1)) +
      
      new_scale_fill() +
      geom_point(aes(y=Incidence_med, group=Sex, fill = Type), shape = 21, size=1, position=position_dodge(1)) +
      geom_linerange(aes(ymin=Incidence_lwr, ymax=Incidence_upr, group=Sex), position=position_dodge(1), linewidth = 0.5) +
      scale_fill_manual(values = fill_vals, guide = "none") +
      # scale_shape_manual(values = shape_vals, guide = NULL) +
      # scale_y_continuous(breaks = seq(0, 50, 10)) +
      scale_x_continuous(breaks=unique(state_df$Age_bin_idx), labels = unique(state_df$Age)) +
      labs (title=state_region_match$State[i], x="Age", y="Incidence (cases per 100,000 pop.)", fill="Sex") +
      theme(axis.text = element_text(colour="black"), axis.ticks = element_line(colour="black"),
            panel.grid.minor.x = element_blank(),
            legend.position = "bottom", legend.position.inside = c(0.19, 0.87),
            legend.background = element_rect(fill = "white", color = "#EBEBEB"),
            axis.text.x = element_text(angle=90, vjust=0.5, hjust=1),
            legend.title = element_blank()) 
    ggsave(here("4_Output_Figures", "Incidence", paste0(state_region_match$State[i],".jpeg")), p, height = 6, width = 10, units = "in")
  }
  
  p_cum_cases <- p_cases %>%
    group_by(Sex, State_code_stan, Age_bin_idx, Type) %>%
    summarise(across(c(med, lwr, upr, Cases), sum)) %>%
    left_join((cum_incidence_state_df %>% select(c(Age_bin_idx, Sex, State_code_stan, Pop))), by=join_by(Age_bin_idx, Sex, State_code_stan)) %>%
    mutate(State = state_region_match$State[match(State_code_stan, state_region_match$State_code_stan)],
           Age = age_bin_match$Age[match(Age_bin_idx, age_bin_match$Age_bin_idx)],
           Incidence = Cases*1e5/11/Pop,
           across(c(med, lwr, upr), ~.x*1e5/11/Pop, .names = "Incidence_{.col}"))
  
  p <- ggplot(p_cum_cases, aes(x=Age_bin_idx)) + 
    facet_wrap(. ~ State, ncol=5, scales="fixed") +
    geom_bar(aes(y = Incidence,  fill=Sex), stat="identity", position="dodge") +
    scale_fill_manual(values=colorpal) +
    theme(legend.position = "bottom", axis.text.x = element_text(angle=90, vjust=0.5, hjust=1, size = 9)) +
    # scale_shape_manual(values = shape_vals, guide = NULL) +
    new_scale_fill() +
    geom_point(aes(y=Incidence_med, group=Sex, fill = Type), shape = 21, size=1, position=position_dodge(1)) +
    geom_linerange(aes(ymin=Incidence_lwr, ymax=Incidence_upr, group=Sex), position=position_dodge(1), linewidth = 0.5) +
    scale_fill_manual(values = fill_vals, guide = "none") +
    scale_x_continuous(breaks=unique(cum_incidence_state_df$Age_bin_idx), labels = unique(cum_incidence_state_df$Age)) +
    labs (title=NULL, x="Age", y="Incidence (cases per 100,000 pop.)", fill="Sex") +
    theme(panel.grid.minor.x = element_blank())
  
  
  ggsave(here("4_Output_Figures", "Supp_Figure_SX_incidence_sex_age_fit_states.pdf"), p, width=8, height = 11, units="in")
  
  p_cum_cases_br <- p_cum_cases %>%
    group_by(Sex, Age_bin_idx, Type) %>%
    summarise(across(c(med, lwr, upr, Cases), sum)) %>%
    left_join((cum_incidence_br_df %>% select(c(Age_bin_idx, Sex, Pop))), by=join_by(Age_bin_idx, Sex)) %>%
    mutate(Age = age_bin_match$Age[match(Age_bin_idx, age_bin_match$Age_bin_idx)],
           Incidence = Cases*1e5/11/Pop,
           across(c(med, lwr, upr), ~.x*1e5/11/Pop, .names = "Incidence_{.col}"))
  
  p <- ggplot(p_cum_cases_br, aes(x=Age_bin_idx)) +
    geom_bar(aes(y=Incidence, fill=Sex), stat="identity", position="dodge") +
    scale_fill_manual(values=colorpal) +
    guides(fill = guide_legend(nrow = 1)) +
    new_scale_fill() +
    geom_point(aes(y=Incidence_med, group=Sex, fill = Type), shape = 21, size=2, position=position_dodge(1)) +
    geom_linerange(aes(ymin=Incidence_lwr, ymax=Incidence_upr, group=Sex), linewidth=0.7, position=position_dodge(1)) +
    scale_fill_manual(values = fill_vals, guide = "none") +
    # scale_shape_manual(values = shape_vals, guide = NULL) +
    scale_y_continuous(breaks = seq(0, 40, 10), limits = c(0, 45)) +
    scale_x_continuous(breaks=unique(cum_incidence_br_df$Age_bin_idx), labels = unique(cum_incidence_br_df$Age)) +
    labs (title=NULL, x="Age", y="Average annual incidence\n(cases per 100,000 pop.)", fill="Sex") +
    theme(axis.text = element_text(colour="black"), axis.ticks = element_line(colour="black"),
          legend.position = "inside", legend.position.inside = c(0.19, 0.87),
          legend.background = element_rect(fill = "white", color = "#EBEBEB"),
          axis.text.x = element_text(angle=90, vjust=0.5, hjust=1),
          panel.grid.minor.x = element_blank(),
          legend.title = element_blank())
  
  ggsave(here("4_Output_Figures", "Fig_1D_incidence_sex_age_with_fit.pdf"), p, width=3.6, height = 2.4, units="in")
  
  p_cases <- p_cases %>%
    filter(Type == "Model")
  max_y <- max(c(p_cases$Cases, p_cases$upr))
  
  p <- ggplot(p_cases, aes(x=Cases, y=med)) +
    geom_point(size=2, shape = 16, alpha = 0.2) +
    geom_errorbar(aes(ymin=lwr, ymax=upr), linewidth=0.7, width=0, alpha = 0.2) +
    geom_line(inherit.aes = F, aes(x=seq(0,max_y,length.out =nrow(p_cases)), y=seq(0,max_y,length.out =nrow(p_cases))), 
              color="red", linetype = "dashed", linewidth = 1) +
    theme(legend.position = "none") +
    scale_x_sqrt() +
    scale_y_sqrt() +
    labs(x="Observed cases", y="Estimated cases")
  
  ggsave(here("4_Output_Figures", "Fig_2B_case_count_fit.pdf"), p, height = 1.8, width=3.25, units = "in")

  p_cases_state <- p_cases %>%
    group_by(State_code_stan, Year) %>%
    summarise(across(.cols=c(med, lwr, upr, Cases, Pop),
                     .fns=sum)) %>%
    mutate(State = state_region_match$State[match(State_code_stan, state_region_match$State_code_stan)]) %>%
    mutate(Incidence = Cases*1e4/Pop,
           across(c(med, lwr, upr), ~.x*1e4/Pop, .names = "Incidence_{.col}"))
  
  p <- p_cases_state %>%
    ggplot(aes(x=Year)) +
    facet_wrap(. ~ State, ncol = 3, scales= "free_y") +
    geom_bar(aes(y = Cases), stat="identity", fill="#4692b0") +
    geom_point(aes(y=med), size=2) +
    geom_linerange(aes(ymin=lwr, ymax=upr), linewidth=0.5) +
    labs(y="Case count") +
    scale_y_continuous(trans = "sqrt") +
    scale_x_continuous(breaks = seq(2014, 2024, by=2)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
  ggsave(here("4_Output_Figures","Supp_Figure_SX_case_fit_year_state.pdf"), p, width=8, height = 11, units="in")
  
  p <- p_cases_state %>%
    ggplot(aes(x=Year)) +
    facet_wrap(. ~ State, ncol = 3, scales= "free_y") +
    geom_bar(aes(y = Incidence), stat="identity", fill="#b75347") +
    geom_point(aes(y=Incidence_med), size=2) +
    geom_linerange(aes(ymin=Incidence_lwr, ymax=Incidence_upr), linewidth=0.5) +
    labs(y="Incidence (cases per 10,000 pop.)") +
    scale_y_continuous(trans = "sqrt") +
    scale_x_continuous(breaks = seq(2014, 2024, by=2)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
  ggsave(here("4_Output_Figures","Supp_Figure_SX_incidence_fit_year_state.pdf"), p, width=8, height = 11, units="in")
  return(NULL)
}

plot_mortality_fit <- function(chik_deaths, age_bin_match){
  pDeaths <- read.csv(here("3_Output_Data", "p_deaths.csv")) %>% 
    mutate(Type = "Model") %>% select(-variable)
  fill_vals <- "black"
  
  if(min(pDeaths$Age_bin_idx) == 2){
    pDeaths <- bind_rows(pDeaths,
                         chik_deaths %>%
                           filter(Age_bin_idx == 1) %>%
                           rename(med = Deaths) %>%
                           select(c(med, Sex, Age_bin_idx, Age)) %>%
                           mutate(Type = "Data")) %>%
      mutate(Type = factor(Type, levels = c("Model", "Data")))
    fill_vals <- c(fill_vals, "white")
  }
  
  pDeaths <- pDeaths %>%
    left_join((chik_deaths %>% select(c(Sex, Age_bin_idx, Deaths, Pop))), by=join_by(Sex, Age_bin_idx)) %>%
    mutate(Death_rate = Deaths*1e5/11/Pop,
           across(c(med, lwr, upr), ~.x*1e5/11/Pop, .names = "Death_rate_{.col}"))
  
  p <- pDeaths %>%
    ggplot(aes(x=Age_bin_idx)) +
    geom_bar(aes(y=Death_rate, fill = Sex), stat="identity", position = "dodge") +
    scale_fill_manual(values=c("#007aa4", "#e84916")) +
    guides(fill = guide_legend(nrow = 1)) +
    new_scale_fill() +
    geom_point(aes(y=Death_rate_med, group = Sex, fill = Type), shape=21, size=2,  position=position_dodge(1)) +
    scale_fill_manual(values = fill_vals, guide = "none") +
    geom_linerange(aes(ymin=Death_rate_lwr, ymax=Death_rate_upr, group=Sex), linewidth=0.7,  position=position_dodge(1)) +
    labs(y="Average annual mortality\n(deaths per 100,000 pop.)", x="Age") +
    scale_x_continuous(breaks=age_bin_match$Age_bin_idx, labels = age_bin_match$Age) +
    theme(legend.position = "inside", legend.position.inside = c(0.19, 0.87),
          legend.background = element_rect(fill = "white", color = "#EBEBEB"),
          axis.text.x = element_text(angle = 90, vjust=0.5, hjust = 1),
          panel.grid.minor.x = element_blank(),
          legend.title = element_blank())
  
  ggsave(here("4_Output_Figures", "Fig_1E_death_incidence_sex_age_with_fit.pdf"), p, width=3.6, height = 2.4, units="in")
}

plot_sero_fit <- function(sero_data_df){
  meanSero <- read.csv(here("3_Output_Data", "seroFit.csv"))
  names(meanSero) <- c("var", "est_med", "est_lwr", "est_upr")
  sero_fit <- cbind(sero_data_df, meanSero)
  max_y <- max(max(c(sero_fit$est_upr, sero_data_df$upr)), 0.8)

  p <- ggplot(sero_fit, aes(x=seroprev, y=est_med)) + 
    geom_point(size=2, shape = 19) +
    geom_errorbar(aes(ymin=est_lwr, ymax=est_upr), linewidth=0.7, width=0) +
    geom_errorbarh(aes(xmin=lwr, xmax=upr), linewidth = 0.7) +
    geom_line(inherit.aes = F, aes(x=seq(0,max_y,length.out =nrow(sero_fit)), y=seq(0,max_y,length.out =nrow(sero_fit))), 
              color="red", linetype = "dashed", linewidth = 1) +
    scale_x_continuous(labels = scales::percent, limits = c(0, max_y)) +
    scale_y_continuous(labels = scales::percent, limits = c(0, max_y)) +
    labs(x="Observed seroprevalence", y="Estimated seroprevalence")
  
  ggsave(here("4_Output_Figures", "Fig_2A_serofit.pdf"), p, height = 1.8, width=3.25, units = "in")
}

get_plot_lambda <- function(param_est_df, model_name, mean_pop_weights, state_region_match, brazil_map, bahia_constraint){
  param_val <- param_est_df[grep("log_lambda", param_est_df$variable), c("variable", "med", "lwr", "upr", "ess_bulk", "ess_tail")] %>%
    mutate(State_code_stan= readr::parse_number(stringr::str_split_i(variable, ",", 1)),
           Year = readr::parse_number(stringr::str_split_i(variable, ",", 2)) + 2013,
           Age_bin_idx = NA,
           med=exp(med),
           lwr=exp(lwr),
           upr=exp(upr),
           State = state_region_match$State[match(State_code_stan, state_region_match$State_code_stan)],
           Model = model_name)
  
  if(bahia_constraint){
    param_val[which(param_val$State_code_stan== 5 & param_val$Year %in% c(2015,2016, 2017)), c("med", "lwr", "upr")] <- param_val[which(param_val$State_code_stan== 5 & param_val$Year == 2014), c("med", "lwr", "upr")]
  }
  
  mean_lambda_state <- param_val %>%
    group_by(State_code_stan) %>%
    summarise(across(.cols=c("med", "lwr", "upr"),
                     .fns=mean)) %>%
    mutate(Year = NA,
           Age_bin_idx = NA,
           State = state_region_match$State[match(State_code_stan, state_region_match$State_code_stan)],
           Model = model_name)
  
  mean_lambda_state <- mean_lambda_state[,intersect(colnames(param_val), colnames(mean_lambda_state))]
  param_val <- param_val[,intersect(colnames(param_val), colnames(mean_lambda_state))]
  param_val <- rbind(param_val, mean_lambda_state)
  
  mean_lambda_br <- as.data.frame(t(apply(mean_lambda_state[,c("med", "lwr", "upr")], 2, weighted.mean, mean_pop_weights)))
  mean_lambda_br$Lambda <- "National mean"
  mean_lambda_br$Model <- model_name
  
  write.table(param_val, here("3_Output_Data", "FOI_state_year_age.csv"), append=T, row.names=F, col.names=F, sep=",")
  write.table(mean_lambda_br, here("3_Output_Data","National_mean_FOI.csv"), append=T, row.names=F, col.names=F, sep=",")
  
  mean_lambda_state$ISO2 <- state_region_match$ISO3_code[match(mean_lambda_state$State_code_stan,
                                                        state_region_match$State_code_stan)]
  
  brazil_map$med <- mean_lambda_state$med[match(brazil_map$ISO2, mean_lambda_state$ISO2)]
  
  max_foi_plot <- (max(brazil_map$med)%/%0.05 + 1)*0.05
  breaks_plot <- seq(0.025, max_foi_plot, 0.025)
  
  p <- ggplot(brazil_map) + 
    geom_sf(aes(fill=med), color = "white", linewidth = 0) +
    theme_void(base_size = 14) +
    scale_fill_distiller(palette="Reds", direction=1, 
                         labels = scales::percent, breaks = breaks_plot, limits = c(0, max_foi_plot)) +
    labs ( x=NULL, y=NULL, fill="Mean FOI") +
    theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(),
          axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(),
          legend.position = "inside", legend.position.inside = c(0.2, 0.25))
  ggsave(here("4_Output_Figures", "Fig_2C_map_mean_foi.pdf"), p, height = 5, width=3.5, units = "in")

  # Rearrange matrix so heatmap is organised geographically, from North to South
  foi_mat <- matrix(param_val$med[!is.na(param_val$Year)], nrow = 27, ncol = 11)
  rownames(foi_mat) <- param_val$State[1:27]
  colnames(foi_mat) <- 2014:2024
  zscores_foi <- t(apply(foi_mat, 1, function(x) (x - mean(x))/sd(x)))
  
  sf_use_s2(FALSE)
  centroids <- st_centroid(brazil_map)
  state_region_match$coords <- centroids$geometry[match(state_region_match$ISO3_code, centroids$ISO2)]
  state_region_match$Lat <- st_coordinates(state_region_match$coords)[,2]
  state_region_match <- state_region_match %>%
    group_by(Region) %>%
    arrange(desc(Lat)) %>%
    ungroup() %>%
    arrange(factor(Region, levels = c("North", "Northeast", "Central-West", "Southeast", "South")))
  
  zscores_foi_df <- as.data.frame(zscores_foi) %>%
    mutate(State = rownames(zscores_foi)) %>%
    tidyr::pivot_longer(cols = where(is.numeric), values_to = "z_score", names_to = "Year")
  
  foi_df <- as.data.frame(foi_mat) %>%
    mutate(State = rownames(foi_mat)) %>%
    tidyr::pivot_longer(cols = where(is.numeric), values_to = "FOI", names_to = "Year")
  
  p1 <- foi_df %>%
    mutate(FOI = if_else(FOI < 0.01, 0.01, FOI)) %>%
    mutate(State = factor(State, levels = rev(state_region_match$State))) %>%
    ggplot(aes(x = Year, y = State, fill = FOI)) +
    geom_tile() +
    scale_fill_distiller(palette="Reds", direction=1, 
                         breaks = c(0.01, 0.03, 0.1, 0.3),
                         labels = c("<1%", "3%", "10%", "30%"), #scales::percent,
                         transform = "log10"
    ) +
    labs(x = NULL, y = NULL, fill = "Annual FOI") +
    scale_x_discrete(breaks = seq(2014, 2024, 2)) +
    theme(legend.position = "bottom",
          axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          axis.text.x = element_text(angle = 90, vjust = 0.5),
          legend.key.width = unit(0.4, "in")) +
    guides(fill = guide_colorbar(title.position="top", title.hjust = 0.5))
  
  p2 <- zscores_foi_df %>%
    mutate(State = factor(State, levels = rev(state_region_match$State))) %>%
    ggplot(aes(x = Year, y = State, fill = z_score)) +
    geom_tile() +
    scale_fill_distiller(palette="BrBG", direction=-1, limits = c(-3.02, 3.02)) +
    # scale_y_discrete(position = "right") +
    scale_x_discrete(breaks = seq(2014, 2024, 2)) +
    labs(x = NULL, y = NULL, fill = "Z-score of annual FOI per state") +
    scale_y_discrete(position = "right") +
    theme(legend.position = "bottom",
          axis.ticks = element_blank(),
          axis.text.x = element_text(angle = 90, vjust = 0.5),
          legend.key.width = unit(0.4, "in")) +
    guides(fill = guide_colorbar(title.position="top", title.hjust = 0.5))
  
  p <- p1 + p2 + plot_annotation(tag_levels = 'A')
  
  ggsave(here("4_Output_Figures", "Supp_Fig_S4_heatmaps_FOI_zscores.png"), p, width = 6, height = 6, units = "in", dpi = 300)
  return(NULL)
}

get_plot_rho_ifr_rr <- function(seed, param_est_df, model_name, case_pop_df, age_bin_match, chik_deaths, dengue_rr_df){
  reporting_rates <-  read.csv(here("3_Output_Data", "mean_rho_sa.csv")) %>% mutate(Type = "Model")
  ifr_df <- param_est_df[grep("log_ifr", param_est_df$variable), c("variable", "med", "lwr", "upr")]
  ifr_df[,c("Sex", "Age_bin_idx")] <- sapply(1:2, function(x) readr::parse_number(stringr::str_split_i(ifr_df$variable, ",", x)))
  
  if(max(ifr_df$Age_bin_idx) == 15) ifr_df$Age_bin_idx <- ifr_df$Age_bin_idx + 1
  ifr_df <- ifr_df %>%
    mutate(Sex = ifelse(Sex == 1, "F", "M"),
           Age = age_bin_match$Age[match(Age_bin_idx, age_bin_match$Age_bin_idx)],
           Type = "Model") %>%
    mutate(across(c(med, lwr, upr), exp)) %>%
    select(-c(variable))
  
  color_pal_rho <- c("#00b9a7", "#ffce00")
  fill_pal_rho <- c("#00b9a7", "#ffce00")
  
  color_pal_ifr <- c("#007aa4", "#e84916")
  fill_pal_ifr <- c("#007aa4", "#e84916")
  
  fill_pal_rr <- c("#FDC300FF", "#0627B3FF", "#BE245BFF")
  col_pal_rr <- c("#FDC300FF", "#0627B3FF", "#BE245BFF")
  
  rr_disease_df <- read.csv(here("3_Output_Data","rr_disease_a.csv")) %>% mutate(Type = "Model")
  rr_death_df <- read.csv(here("3_Output_Data","rr_death_a.csv")) %>% mutate(Type = "Model")
  
  if(min(reporting_rates$Age_bin_idx) == 2){
    log_lambda_draws_df <- read.csv(here("3_Output_Data","log_lambda_draws.csv"))
    
    set.seed(seed)
    sampled_df <- data.frame(.draw = sample(max(log_lambda_draws_df$.draw), 1e3, replace = T))
    
    log_lambda_draws_df <- log_lambda_draws_df %>%
      right_join(sampled_df, by=".draw", relationship = "many-to-many") %>%
      mutate(foi = exp(value)) %>%
      select(-c(.chain, .iteration, variable, value))
    
    rg_df <- read.csv(here("3_Output_Data","Rg_draws.csv")) %>%
      right_join(sampled_df, by=".draw", relationship = "many-to-many") %>%
      filter(Sex == "F", Age_bin_idx %in% 5:10) %>%
      left_join((case_pop_df %>% 
                   select(c(Sex, State_code_stan, Year, Age_bin_idx, Pop))),
                by = join_by(Sex, State_code_stan, Year, Age_bin_idx)) %>%
      group_by(State_code_stan, Year, .draw) %>%
      summarise(Imm_mat = weighted.mean(value, Pop)) %>%
      mutate(Year = Year + 1,
             Imm_mat = ifelse(Year == 2025, 0, Imm_mat),
             Year = ifelse(Year == 2025, 2014, Year)) %>% 
      arrange(State_code_stan, Year, .draw)
    
    infections_infants_df <- log_lambda_draws_df %>%
      left_join(rg_df, by = join_by(State_code_stan, Year, .draw))  %>%
      slice(rep(1:n(), each = 2)) %>%
      mutate(Sex = rep(c("F", "M"), times = nrow(log_lambda_draws_df))) %>%
      left_join((case_pop_df %>% 
                   filter(Age_bin_idx == 1) %>%
                   select(c(Sex, State_code_stan, Year, Cases, Pop))), by = join_by(Sex, State_code_stan, Year)) %>%
      mutate(Infections = ((1-Imm_mat) + Imm_mat*0.5)*foi*Pop) 
    
    mean_rho_v1_df <- read.csv(here("3_Output_Data","mean_rho.csv"))
    tot_infections_v1_df <- read.csv(here("3_Output_Data","total_infections.csv")) %>%
      filter(Age_bin_idx >= 2) %>%
      summarise(across(c(med,lwr,upr), sum, .names = "Inf_{.col}"))
    
    mean_rho_v1_df <- bind_cols(mean_rho_v1_df,
                                tot_infections_v1_df) %>%
      mutate(Group = "Older")
    
    tot_cases_infants <- case_pop_df %>%
      filter(Age_bin_idx == 1) %>%
      summarise(Cases = sum(Cases))
    
    mean_rho_infants_df <- infections_infants_df %>%
      group_by(.draw) %>%
      summarise(Infections = sum(Infections)) %>%
      ungroup() %>%
      summarise(Inf_med = median(Infections),
                Inf_lwr = quantile(Infections, prob = 0.025),
                Inf_upr = quantile(Infections, prob = 0.975),
                med = tot_cases_infants$Cases/Inf_med,
                lwr = tot_cases_infants$Cases/Inf_upr,
                upr = tot_cases_infants$Cases/Inf_lwr,
                Group = "Infants")
    
    mean_rho_overall <- bind_rows(mean_rho_v1_df,
                                  mean_rho_infants_df) %>%
      summarise(med = weighted.mean(med, Inf_med),
                lwr = weighted.mean(lwr, Inf_upr),
                upr = weighted.mean(upr, Inf_lwr))
    
    write.csv(mean_rho_overall, here("3_Output_Data","mean_rho_overall.csv"))
    
    rep_infants_df <- infections_infants_df %>%
      group_by(.draw, Sex) %>%
      summarise(rho = sum(Cases)/sum(Infections)) %>%
      tidyr::pivot_wider(names_from = Sex, values_from = rho, names_prefix = "rho_") %>%
      ungroup() %>%
      mutate(rel_risk = rho_M/rho_F) %>%
      reframe(across(c(rho_F, rho_M, rel_risk), ~ quantile(.x, c(0.5, 0.025, 0.975))))
    
    rep_infants_df <- as.data.frame(t(rep_infants_df))
    names(rep_infants_df) <- c("med", "lwr", "upr")
    rownames(rep_infants_df) <- NULL
    rep_infants_df <- rep_infants_df %>%
      mutate(Sex = c("F", "M", NA),
             Age_bin_idx = 1, 
             Age = "<1",
             Type = "Post-hoc")
    
    reporting_rates <- bind_rows(reporting_rates,
                                 rep_infants_df %>%
                                   filter(!is.na(Sex))) %>%
      mutate(Type = factor(Type, levels = c("Model", "Post-hoc")))
    
    rr_disease_df <- bind_rows(rr_disease_df,
                               rep_infants_df %>%
                                 filter(is.na(Sex)) %>%
                                 select(-Sex)) %>%
      mutate(Type = factor(Type, levels = c("Model", "Post-hoc")))
    
    mean_ifr_v1_df <- read.csv(here("3_Output_Data","mean_ifr.csv"))
    mean_ifr_v1_df <- bind_cols(mean_ifr_v1_df,
                                tot_infections_v1_df) %>%
      mutate(Group = "Older")
    
    tot_deaths_infants <- chik_deaths %>%
      filter(Age_bin_idx == 1) %>%
      summarise(Deaths = sum(Deaths))
    
    mean_ifr_infants_df <- infections_infants_df %>%
      group_by(.draw) %>%
      summarise(Infections = sum(Infections)) %>%
      ungroup() %>%
      summarise(Inf_med = median(Infections),
                Inf_lwr = quantile(Infections, prob = 0.025),
                Inf_upr = quantile(Infections, prob = 0.975),
                med = tot_deaths_infants$Deaths/Inf_med,
                lwr = tot_deaths_infants$Deaths/Inf_upr,
                upr = tot_deaths_infants$Deaths/Inf_lwr,
                Group = "Infants")
    
    mean_ifr_overall <- bind_rows(mean_ifr_v1_df,
                                  mean_ifr_infants_df) %>%
      summarise(med = weighted.mean(med, Inf_med),
                lwr = weighted.mean(lwr, Inf_upr),
                upr = weighted.mean(upr, Inf_lwr))
    
    write.csv(mean_ifr_overall, here("3_Output_Data", "mean_ifr_overall.csv"))
    
    ifr_infants_df <- infections_infants_df %>%
      group_by(.draw, Sex) %>%
      summarise(Infections = sum(Infections)) %>%
      left_join(chik_deaths %>%
                  filter(Age_bin_idx == 1) %>%
                  select(Sex, Deaths), by = "Sex") %>%
      mutate(ifr = Deaths/Infections) %>%
      tidyr::pivot_wider(id_cols = .draw, names_from = Sex, values_from = ifr, names_prefix = "ifr_") %>%
      ungroup() %>%
      mutate(rel_risk = ifr_M/ifr_F) %>%
      reframe(across(c(ifr_F, ifr_M, rel_risk), ~ quantile(.x, c(0.5, 0.025, 0.975))))
    
    ifr_infants_df <- as.data.frame(t(ifr_infants_df))
    names(ifr_infants_df) <- c("med", "lwr", "upr")
    rownames(ifr_infants_df) <- NULL
    ifr_infants_df <- ifr_infants_df %>%
      mutate(Sex = c("F", "M", NA),
             Age_bin_idx = 1, 
             Age = "<1",
             Type = "Post-hoc")
    
    ifr_df <- bind_rows(ifr_df,
                        ifr_infants_df %>%
                          filter(!is.na(Sex))) %>%
      mutate(Type = factor(Type, levels = c("Model", "Post-hoc")))
    
    rr_death_df <- bind_rows(rr_death_df, 
                             ifr_infants_df %>%
                               filter(is.na(Sex)) %>%
                               select(-Sex)) %>%
      mutate(Type = factor(Type, levels = c("Model", "Post-hoc")))
    
    fill_pal_rho <- c(fill_pal_rho, "white", "white")
    fill_pal_ifr <- c(fill_pal_ifr, "white", "white")
    fill_pal_rr <- c(fill_pal_rr, "white", "white", "white")
    
  }
  
  p <- reporting_rates %>%
    ggplot(aes(x=Age_bin_idx, y=med, color=Sex)) +
    geom_linerange(aes(ymin=lwr, ymax=upr), position = position_dodge2(width=0.65), linewidth = 0.7) +
    geom_point(aes(fill = interaction(Sex,Type)), shape = 21, size=2, position=position_dodge2(width = 0.65)) +
    scale_colour_manual(values=color_pal_rho) +
    scale_y_continuous(trans="log10", 
                       # breaks = c(0.003, 0.007, 0.01, 0.013, 0.017), 
                       labels = scales::percent, 
                       limits = c(min(reporting_rates$lwr) - 0.001, max(reporting_rates$upr) + 0.001)) +
    scale_x_continuous(breaks=unique(reporting_rates$Age_bin_idx), labels=unique(reporting_rates$Age)) +
    scale_fill_manual(values = fill_pal_rho, guide = NULL) +
    labs(y="Marginal probability of detection", x="Age") +
    theme(axis.text.x = element_text(angle=90, hjust=1, vjust = 0.5), 
          legend.position = "inside",
          legend.position.inside = c(0.75, 0.15),
          legend.background = element_rect(fill = "white", color = "#EBEBEB"),
          panel.grid.minor.x = element_blank(),
          legend.title = element_blank()) +
    guides(color = guide_legend(nrow = 1, override.aes = list(shape = 19)))
  ggsave(here("4_Output_Figures", "Fig_3B_prob_reporting.pdf"), p, width=3.25, height = 3, units="in")
  write.csv(reporting_rates, here("3_Output_Data", "mean_rho_sa_vf.csv"), row.names = F)
  
  p <- ifr_df %>%
    ggplot(aes(x=Age_bin_idx, y=med, color=Sex)) +
    geom_linerange(aes(ymin=lwr, ymax=upr), position = position_dodge2(width=0.65), linewidth = 0.7) +
    geom_point(aes(fill = interaction(Sex, Type)), shape = 21, size=2, position=position_dodge2(width = 0.65)) +
    scale_colour_manual(values=color_pal_ifr) +
    scale_fill_manual(values = fill_pal_ifr, guide = NULL) +
    scale_y_continuous(trans="log10", labels = scales::percent, 
                       breaks = c(3e-6, 1e-5, 3e-5,  1e-4, 3e-4, 5e-4)) +
    scale_x_continuous(breaks=unique(ifr_df$Age_bin_idx), labels=unique(ifr_df$Age)) +
    labs(y="Infection Fatality Rate", x="Age") +
    theme(axis.text.x = element_text(angle=90, hjust=1, vjust = 0.5), 
          legend.position="inside", legend.position.inside = c(0.75, 0.15),
          panel.grid.minor.x = element_blank(),
          legend.title = element_blank(),
          legend.background = element_rect(fill = "white", color = "#EBEBEB")) +
    guides(color = guide_legend(nrow = 1, override.aes = list(shape = 19)))
  ggsave(here("4_Output_Figures","Fig_3C_ifr_sex_age.pdf"), p, width=3.25, height = 3, units="in") 
  write.csv(ifr_df, here("3_Output_Data",  "mean_ifr_sa.csv"), row.names = F)
  
  rr_tot_df <- bind_rows(rr_death_df %>% mutate(Info = "Mortality"), 
                         rr_disease_df %>% mutate(Info = "Detection")) %>%
    mutate(across(c(med,lwr,upr), ~1/.x))
  
  chik_adj_rr_df <- rr_disease_df %>%
    left_join(dengue_rr_df %>% select(Age, RR) %>% rename(RR_dengue = RR), by = "Age") %>%
    mutate(across(c(med, lwr, upr), ~ (1/.x)/RR_dengue)) %>%
    mutate(Info= "Detection - adjusted") %>%
    select(Age, Age_bin_idx, med, lwr, upr, Info, Type)
  
  rr_plot_df <- bind_rows(rr_tot_df,
                          chik_adj_rr_df)
  
  p <- rr_plot_df %>%
    mutate(Info = factor(Info, levels = c("Detection", "Detection - adjusted",  "Mortality"))) %>%
    ggplot(aes(x=Age_bin_idx, y=med, color = Info, group = Info)) +
    geom_hline(yintercept = 1, color="red", linetype = "dashed", linewidth = 1) +
    geom_point(aes(fill = interaction(Info, Type)), shape = 21, size=1.5, position = position_dodge2(width = 0.5)) +
    geom_linerange(aes(ymin=lwr, ymax=upr), linewidth = 0.5, position = position_dodge2(width = 0.5)) +
    geom_line(linewidth = 0.5, position = position_dodge2(width = 0.5)) +
    scale_y_continuous(trans="log10", breaks = c(0.3, 1.0, 3.0, 5.0), limits = c(0.22, 5.29)) +
    scale_x_continuous(breaks=unique(rr_tot_df$Age_bin_idx), labels=unique(rr_tot_df$Age)) +
    scale_fill_manual(values = fill_pal_rr, guide = NULL) +
    scale_color_manual(values = col_pal_rr) +
    labs(y="RR (F : M)", x="Age", color = NULL, fill = NULL) +
    theme(axis.text.x = element_text(angle=90, hjust=1, vjust = 0.5),
          panel.grid.minor.x = element_blank(),
          legend.text = element_text(size = txt_size-2),
          legend.position="bottom") +
          # legend.position.inside = c(0.85, 0.9),
          # legend.background = element_rect(fill = "white", color = "grey90")) +
    guides(color = guide_legend(override.aes = list(shape = 19), nrow = 1))
  ggsave(here("4_Output_Figures", "Fig_3D_rr_detection_death_FvsM.pdf"), p, height = 3.5, width = 3, units = "in")
  
  return(NULL)
}

get_rho_t <- function(){
  reporting_rates <-  read.csv(here("3_Output_Data", "mean_rho_t.csv"))
  
  p <- reporting_rates %>%
    ggplot(aes(x=Year, y=med)) +
    geom_linerange(aes(ymin=lwr, ymax=upr), linewidth = 0.7, color ="#ca2360") +
    geom_point(size = 2, color ="#ca2360") +
    scale_y_continuous(labels = scales::percent, transform = "log10") +
    scale_x_continuous(breaks=2014:2024) +
    labs(y="Marginal probability of detection", x="Year") +
    theme(axis.text.x = element_text(angle=90, hjust=1, vjust = 0.5),
          panel.grid.minor.x = element_blank())
  ggsave(here("4_Output_Figures", "Fig_3A_prob_reporting_time.pdf"), p, width=3.25, height = 3, units="in")
  
  return(NULL)
}

get_prop_infected <- function(case_pop_df, state_region_match, brazil_map){
  top_states <- case_pop_df %>%
    group_by(State_code_stan) %>%
    summarise(tot_cases = sum(Cases)) %>%
    arrange(desc(tot_cases)) %>% 
    slice_head(n=10) %>% select(State_code_stan) %>% as.vector()
  
  rg_df <- read.csv(here("3_Output_Data", "Rg.csv")) %>%
    left_join((case_pop_df %>% select(Sex, State_code_stan, Age_bin_idx, Year, Pop)), by = join_by(Sex, State_code_stan, Age_bin_idx, Year)) %>%
    mutate(Top_10 = State_code_stan %in% top_states[[1]],
           Top_5 = State_code_stan %in% top_states[[1]][1:5],
           Region = state_region_match$Region[match(State_code_stan, state_region_match$State_code_stan)])
  
  rg_top_5_df <- rg_df %>%
    filter(Top_5) %>%
    group_by(Year) %>%
    summarise(across(c(med, lwr, upr), ~ weighted.mean(.x, Pop))) %>%
    mutate(Type = "Top 5 States")
  
  rg_all_df <- rg_df %>%
    group_by(Year) %>%
    summarise(across(c(med, lwr, upr), ~ weighted.mean(.x, Pop))) %>%
    mutate(Type = "Brazil")
  
  rg_plot_df <- bind_rows(rg_top_5_df, rg_all_df)

  colpal <- c( "#C8350DFF",  "#F9C000FF")
  
  p <- rg_plot_df %>%
    mutate(Type = factor(Type, levels = c("Top 5 States", "Top 10 States", "Brazil"))) %>%
    ggplot(aes(x = Year, y = med)) +
    geom_ribbon(aes(fill = Type, ymin = lwr, ymax = upr), alpha = 0.6) +
    geom_line(aes(color = Type), linewidth = 0.7) +
    scale_color_manual(values = colpal) +
    scale_fill_manual(values = colpal) +
    scale_y_continuous(breaks =seq(0, 0.35, 0.05), limits = c(0, 0.35), labels = scales::percent) +
    scale_x_continuous(breaks = seq(2014, 2024, 2), labels = seq(2014, 2024, 2)) +
    labs(fill = NULL, color = NULL, y= "Proportion Infected") +
    theme(legend.position = "inside",
          legend.position.inside = c(0.85, 0.15),
          panel.grid.minor.y = element_blank(),
          legend.background = element_rect(fill = "white", color = "#EBEBEB"))
  ggsave(here("4_Output_Figures", "Fig_2E_infection_time.pdf"), p, width=6, height = 3, units="in")
  
  rg_2024_df <- rg_df %>%
    filter(Year == 2024) %>%
    group_by(State_code_stan) %>%
    summarise(across(c(med, lwr, upr), ~ weighted.mean(.x, Pop))) %>%
    mutate(ISO2 = state_region_match$ISO3_code[match(State_code_stan, state_region_match$State_code_stan)])
  
  brazil_map <-  brazil_map %>%
    left_join((state_region_match %>% select(ISO3_code, State, State_code_stan)), by = c("ISO2"="ISO3_code")) %>%
    mutate(Immunity = rg_2024_df$med[match(ISO2, rg_2024_df$ISO2)],
           Top_5 = State_code_stan %in% top_states[[1]][1:5]) %>% arrange(Top_5)
  
  max_imm_plot <- (max(brazil_map$Immunity)%/%0.05 + 1)*0.05
  
  if(max_imm_plot == 0.8){
    breaks_plot <- seq(0, max_imm_plot, length.out = 5)
  } else if(max_imm_plot == 0.75){
    breaks_plot <- seq(0, max_imm_plot, length.out = 4)
  } else {
    breaks_plot <- round(seq(0, max_imm_plot, length.out = 4), 2)
  } 
  
  p <- ggplot(brazil_map) + 
    geom_sf(aes(fill=Immunity, color = Top_5, linewidth = Top_5)) +
    theme_void() +
    scale_color_manual(values = c("#C8350DFF", "#333333"),
                       labels = c("Top 5 states"),
                       breaks = c(TRUE), na.value = "white") +
    scale_linewidth_manual(values = c(0, 0.5), na.value = 0, guide = NULL) +
    scale_fill_distiller(palette="BuGn", direction=1, labels = scales::percent,
                         breaks = breaks_plot, limits = c(0, max_imm_plot)) +
    labs ( x=NULL, y=NULL, fill="Proportion\ninfected", linewidth = NULL, color = NULL) +
    theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(),
          axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(),
          legend.position = "inside", legend.position.inside = c(0.2, 0.25))
  ggsave(here("4_Output_Figures", "Fig_2D_map_prop_infected_2024.pdf"), p,  height = 5, width=3.5, units = "in")
}
