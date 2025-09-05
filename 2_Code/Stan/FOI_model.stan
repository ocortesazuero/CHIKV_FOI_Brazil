//----- Time-varying FOI model -----//
  data {
    int nA; // N age groups
    int nT; // N time points
    int nL; // N locations
    int nSero; // N serostudies
    int nCasePoints; // N case points - used in log likelihood output
    int nPoints; // N total points - used in log likelihood output
    array[2,nL,nT,nA] int cases; // Cases per sex, location, time point, age group
    array[2,nL] matrix[nT,nA] pop; // Population per cohort
    array[2,nL,nT,nA] int pointIndex; // Indexing for log likelihood output
    array[nA] int aMin; // Index for age groups
    array[nA] int aMax; // Index for age groups
    array[nSero] int seroYear; // Characteristics of serostudies
    array[nSero] int seroLocs;
    array[nSero] int seroAgeMin;
    array[nSero] int seroAgeMax;
    array[nSero] int seroPos;
    array[nSero] int seroSize;
    array[nSero] int seroSex;
    array[2,nA] int deaths; // Deaths per sex, age group
    array[2,nA] int pointDeathIndex; // Indexing for log likelihood output
  }

parameters {
  array[nL,nT] real<upper=-0.1> log_lambda; // Log FOI per location, year
  array[2,nA] real logit_rho; // Sex- and age-dependent reporting probabilities
  ordered[nT] logit_rho_t; // Time-dependent reporting probabilities
  real<lower=1> phi; // Overdispersion parameter for negative binomial distribution
  array[2,nA] real<upper=-2> log_ifr; // Log IFR per sex and age group
}

transformed parameters {
  array[2,nL] matrix<lower=0, upper=1>[nT,100] S; // Susceptible per sex, loc, time, individual age
  array[2,nL] matrix<lower=0, upper=1>[nT,100] I; // Infected
  array[2,nL] matrix<lower=0, upper=1>[nT,100] R; // Immune
  
  array[2,nL] matrix<lower=0, upper=1>[nT,nA] Sg; // Susceptible per sex, loc, time, aggregated age group
  array[2,nL] matrix<lower=0, upper=1>[nT,nA] Ig; // Infected
  array[2,nL] matrix<lower=0, upper=1>[nT,nA] Rg; // Immune
  
  array[2,nL,nT,nA] real pLabCases; // Predicted cases per sex, loc, time, age group
  array[2,nA] real pDeaths; // Predicted deaths per sex, age gorup 
  
  // Re-scaling parameters
  array[nL,nT] real lambda = exp(log_lambda);
  array[2,nA] real rho= inv_logit(logit_rho);
  array[2,nA] real ifr = exp(log_ifr);
  vector<lower = 0, upper=1>[nT] rho_t = inv_logit(logit_rho_t);
  
  array[2,nL,nT,nA] real ITot; // Total number of infections
  array[2,nA] real ITotSexAge = rep_array(0, 2, nA);
  array[nSero] real meanSero; // Expected seroprevalence for each serostudy
  

  for(t in 2:4){
    lambda[5,t] = lambda[5,1]; // Constraint on Bahia FOIs to aid model convergence
  }
  
  // First year of transmission starts with wholly susceptible population
  for(s in 1:2){
    for(l in 1:nL){
      S[s,l,1,] = 1 - lambda[l,1]*rep_row_vector(1,100);
      I[s,l,1,] = lambda[l,1]*rep_row_vector(1,100);
      R[s,l,1,] = lambda[l,1]*rep_row_vector(1,100);
    }
  }
  
  // Loop through subsequent yearly timesteps
  for(s in 1:2){
    for(l in 1:nL){
      for(t in 2:nT){
        S[s,l,t,1] = 1 - lambda[l,t]; // new babies
        I[s,l,t,1] = lambda[l,t];
        
        S[s,l,t,2:100] = S[s,l,t-1,1:99] - lambda[l,t]*S[s,l,t-1,1:99];
        I[s,l,t,2:100] = lambda[l,t]*S[s,l,t-1,1:99];
        R[s,l,t,1:100] = 1-S[s,l,t,1:100];
      }
    }
  }
  
  // Aggregate infections to age groups and report cases
  for(s in 1:2) for(l in 1:nL) for(a in 1:nA) for(t in 1:nT) {
    Sg[s,l,t,a] = mean(S[s,l,t,aMin[a]:aMax[a]]);
    Ig[s,l,t,a] = mean(I[s,l,t,aMin[a]:aMax[a]]);
    Rg[s,l,t,a] = mean(R[s,l,t,aMin[a]:aMax[a]]);
    
    ITot[s,l,t,a] = Ig[s,l,t,a]*pop[s,l,t,a];
    pLabCases[s,l,t,a] = rho[s,a]*rho_t[t]*ITot[s,l,t,a];
    
    ITotSexAge[s,a] += ITot[s,l,t,a];
  }
  
  // Deaths
  for(s in 1:2) for(a in 1:nA) pDeaths[s,a] = ifr[s,a]*ITotSexAge[s,a];
  
  // Seroprevalence in the setting of each serostudy
  for (i in 1:nSero){
    meanSero[i] = sum(Rg[1,seroLocs[i],seroYear[i],seroAgeMin[i]:seroAgeMax[i]].*pop[1,seroLocs[i],seroYear[i],seroAgeMin[i]:seroAgeMax[i]] + seroSex[i]*Rg[2,seroLocs[i],seroYear[i],seroAgeMin[i]:seroAgeMax[i]].*pop[2,seroLocs[i],seroYear[i],seroAgeMin[i]:seroAgeMax[i]])/sum(pop[1,seroLocs[i],seroYear[i],seroAgeMin[i]:seroAgeMax[i]] + seroSex[i]*pop[2,seroLocs[i],seroYear[i],seroAgeMin[i]:seroAgeMax[i]]);
  }
  
}


model {
  
  //--- Priors ---//
    for(l in 1:nL){
      log_lambda[l,] ~ normal(-5,1);
    }
  
  phi ~ normal(4,1);
  for(s in 1:2){
    log_ifr[s,] ~ normal(-8, 2);
    logit_rho[s,] ~ normal(-4,2);
  }
  logit_rho_t ~ normal(-4, 2);
  
  //--- Likelihood ---//
  for (s in 1:2) for (l in 1:nL) for(t in 1:nT) cases[s,l,t,] ~ neg_binomial_2(pLabCases[s,l,t,], phi);

  for(i in 1:nSero){
    seroPos[i] ~ binomial(seroSize[i], meanSero[i]);
  }
  
  for (s in 1:2) deaths[s,] ~ poisson(pDeaths[s,]);
}

generated quantities {
  // Exporting point log likelihood as a one-dimensional array
  // Calculating marginal probabilities of detection according to sex and age, and according to time
  // Calculating risk ratios (M:F) of detection of disease and of death
  
  
  array[nPoints] real log_lik;
  real log_lik_sum = 0.0;
  real rmse = 0.0;
  
  array[2,nA] real mean_rho_sa;
  array[2] real mean_rho_s;
  array[nA] real mean_rho_a;
  array[nT] real mean_rho_t;
  real mean_rho;
  
  array[2] real mean_ifr_s;
  array[nA] real mean_ifr_a;
  real mean_ifr;
  
  array[nA] real rr_disease_a;
  array[nA] real rr_death_a;
  
  array[2,nL,nA] real p_deaths_sla;
  
  real pLabCasesTot = 0;
  real pInfTot = 0;
  real pDeathsTot = 0;
  
  for (s in 1:2){
    for(l in 1:nL){
      for(t in 1:nT){
        for(a in 1:nA){
          log_lik[pointIndex[s,l,t,a]] = neg_binomial_2_lpmf(cases[s,l,t,a] | pLabCases[s,l,t,a], phi);
          log_lik_sum += log_lik[pointIndex[s,l,t,a]];
          rmse += (pLabCases[s,l,t,a] - cases[s,l,t,a])^2;
        }
      }
    }
  }
  
  for (i in 1:nSero){
    log_lik[i + nCasePoints] = binomial_lpmf(seroPos[i] | seroSize[i], meanSero[i]);
    log_lik_sum += log_lik[i + nCasePoints];
    rmse += ((1.0*seroPos[i]/1.0*seroSize[i]) - meanSero[i])^2;
  }
  
  for(s in 1:2){
    for(a in 1:nA){
      log_lik[pointDeathIndex[s,a] + nCasePoints + nSero] = poisson_lpmf(deaths[s,a] | pDeaths[s,a]);
      log_lik_sum += log_lik[pointDeathIndex[s,a] + nCasePoints + nSero];
      rmse += (pDeaths[s,a] - deaths[s,a])^2;
    }
  }
  
  rmse = (rmse/nPoints) ^ .5;
  
  for(a in 1:nA){
    rr_disease_a[a] = rho[2,a]/rho[1,a];
    rr_death_a[a] = ifr[2,a]/ifr[1,a];
  }
  
  for(s in 1:2) for(a in 1:nA){
    real tmpSumCases = 0;
    real tmpSumInf = 0;
    for (t in 1:nT){ 
      tmpSumCases += sum(pLabCases[s,,t,a]);
      tmpSumInf += sum(ITot[s,,t,a]);
    }
    pLabCasesTot += tmpSumCases;
    pInfTot += tmpSumInf;
    mean_rho_sa[s,a] = tmpSumCases/tmpSumInf;
  }
  
  for(s in 1:2){
    real tmpSumCases = 0;
    real tmpSumInf = 0;
    
    for(t in 1:nT) for(a in 1:nA){
      tmpSumCases += sum(pLabCases[s,,t,a]);
      tmpSumInf += sum(ITot[s,,t,a]);
    }
    pDeathsTot += sum(pDeaths[s,]);
    mean_rho_s[s] = tmpSumCases/tmpSumInf;
    mean_ifr_s[s] = sum(pDeaths[s,])/tmpSumInf;
  }
  
  for(a in 1:nA){
    real tmpSumCases = 0;
    real tmpSumInf = 0;
    
    for(s in 1:2) for(t in 1:nT){
      tmpSumCases += sum(pLabCases[s,,t,a]);
      tmpSumInf += sum(ITot[s,,t,a]);
    }
    
    mean_rho_a[a] = tmpSumCases/tmpSumInf;
    mean_ifr_a[a] = sum(pDeaths[,a])/tmpSumInf;
  }
  
  for(t in 1:nT){
    real tmpSumCases = 0;
    real tmpSumInf = 0;
    
    for(s in 1:2) for(a in 1:nA){
      tmpSumCases += sum(pLabCases[s,,t,a]);
      tmpSumInf += sum(ITot[s,,t,a]);
    }
    
    mean_rho_t[t] = tmpSumCases/tmpSumInf;
  }
  
  mean_rho = pLabCasesTot/pInfTot;
  mean_ifr = pDeathsTot/pInfTot;
  
  for(s in 1:2) for(l in 1:nL) for(a in 1:nA){
    p_deaths_sla[s,l,a] = ifr[s,a]*sum(ITot[s,l,,a]);
  }
}


