# R version 4.3.1 (2023-06-16)
# Platform: aarch64-apple-darwin20 (64-bit)
# Running under: macOS Sonoma 14.6.1

# Packages
# tidyverse_2.0.0
# rstan_2.33.1.9000 
# coda_0.19-4
# loo_2.6.0

rm(list=ls())

#### > 1. Initialization ####
#### >> 1.1. Packages ####
libs <- c(
  'tidyverse',
  'rstan', 'bayesplot', 
  'coda', 'loo' 
)
invisible(lapply(libs, library, character.only = T))

#### >> 1.2. Datasets ####
df_hetero_simul <- read_csv("df-hetero-simul.csv")

#### > 2. Model ####
#### >> 2.1. Define model ####
stan_code = '
data{
  int n; // observations
  int m; // species number
  vector[n] decompC; // carbon decomposition
  vector[n] litterBM; // total litter biomass
  vector[n] litterSR; // litter species richness
  matrix[n, m] p; // proportion of species in litter biomass
}

parameters{
  vector[m] b; // species specific proportion parameter
  vector[m] a; // species specific interraction parameter
  real c; // litter biomass parameter
  real d; // litter species richness parameter
  real<lower=0> sigma; // sigma
}  
  
model{
  // auxiliary variables
  real speff;
  real intereff;

  // Priors
  b ~ normal(0,10); // species specific proportion parameter
  a ~ normal(0,10); // species specific interraction parameter
  c ~ normal(0,10); // litter biomass parameter
  d ~ normal(0,10); // litter species richness parameter
  sigma ~ normal(0,10); // sigma
  
  // Likelihood
  
  for(k in 1:n){
    speff = 0.0;
    intereff = 0.0;
    for(i in 1:m){
      // summed species effects
      speff = speff + b[i] * p[k,i];
      // summed interaction effects: use trick for avoiding double sums
      intereff = intereff + a[i]*p[k,i]*(1-p[k,i]);
     }
    decompC[k] ~ normal(speff + intereff + c * litterBM[k] + d * litterSR[k], sigma); 
  }
}
'
stan_model <- stan_model(model_code = stan_code)

#### >> 2.2. Fit model ####
#### >>> 2.2.1 Carbon ####
dat_stan_C.Ma <- list(
  n = nrow(df_hetero_simul),
  m = 16,
  decompC = df_hetero_simul$C.loss_Ma1,
  litterBM = rowSums(df_hetero_simul[,8:23]),
  litterSR = rowSums(df_hetero_simul[,8:23]!=0),
  p = df_hetero_simul[,c(8,14,15,10,9,13,18,16,11,12,17,19,20,21,22,23)]/rowSums(df_hetero_simul[,8:24])
)

fit_carbon <- sampling(stan_model, 
                       data = dat_stan_C.Ma,
                       iter = 3000,
                       warmup = 1000,
                       chains = 3,
                       cores = 3,
                       refresh=10,
                       control=list(max_treedepth=10)
)

#### >>> 2.2.2 Nitrogen ####
dat_stan_N.Ma <- list(
  n = nrow(df_hetero_simul),
  m = 16,
  decompC = df_hetero_simul$N.loss_Ma1,
  litterBM = rowSums(df_hetero_simul[,8:23]),
  litterSR = rowSums(df_hetero_simul[,8:23]!=0),
  p = df_hetero_simul[,c(8,14,15,10,9,13,18,16,11,12,17,19,20,21,22,23)]/rowSums(df_hetero_simul[,8:24]) # each row = 1 observations 
)

fit_nitrogen <- sampling(stan_model, 
                         data = dat_stan_N.Ma,
                         iter = 3000,
                         warmup = 1000,
                         chains = 3,
                         cores = 3,
                         refresh=10,
                         control=list(max_treedepth=10)
)

#### > 3. Predict response from model ####
#### >> 3.1 Function #### 
pred.decomp = function(newdata, fit){
  par = summary(fit)$summary
  pred = 
    sum(par[1:16,1] * newdata[1:16]) + # identity effect
    sum(par[17:32,1] * newdata[1:16] * (1-newdata[1:16])) + # interaction effect
    par[33,1] * newdata[17] + # biomass effect
    par[34,1] * newdata[18] # species richness effect
  return(pred)
}

#### >> 3.2. Fit quality ####
#### >>> 3.2.1. Carbon model ####
#### >>>> 3.2.1.1 Predict original data ####
df = bind_cols(dat_stan_C.Ma[['p']], 
               dat_stan_C.Ma[['litterBM']],
               dat_stan_C.Ma[['litterSR']], 
               dat_stan_C.Ma[['decompC']])
df[,20] = NA
df = as.matrix(df)
for(i in 1: nrow(df)){
  df[i,20] = pred.decomp(newdata = df[i, 1:18], fit = fit_carbon)
}

#### >>>> 3.2.1.2. Bayesian R2 ####
R2.C = var(df[,20])/(var(df[,20]) + var(df[,20] - df[,19]))
R2.C

#### >>>> 3.2.1.3 Predicted vs measured ####
plot(df[,20]~df[,19])
abline(a = 0, b= 1)

#### >>>> 3.2.1.4 Residuals distribution ####
plot(sqrt(abs(df[,20]-df[,19]))~df[,19])

#### >>>> 3.2.1.5 Posterior distribution ####
ggpubr::ggarrange(
  as.matrix(fit_carbon) |>
    data.frame() |>
    select(1:16) |>
    pivot_longer(cols = 1:16) |>
    ggplot(aes(x = value, fill = name)) + 
    geom_density() + 
    geom_vline(xintercept = 0) + 
    facet_grid(rows = vars(name)) + 
    labs(title = 'Species-specific decomposition') + 
    theme_bw() +
    theme(axis.text.y = element_blank(), 
          legend.position = 'none',
          strip.text.y = element_blank()),
  
  as.matrix(fit_carbon) %>%
    data.frame() |>
    select(17:32) |>
    pivot_longer(cols = 1:16) %>%
    ggplot(aes(x = value, fill = name)) + 
    geom_density() + 
    geom_vline(xintercept = 0) + 
    facet_grid(rows = vars(name)) + 
    theme_bw()+
    labs(title = 'Interaction effect') + 
    theme(axis.text.y = element_blank(), 
          legend.position = 'none',
          strip.text.y = element_blank()),
  as.matrix(fit_carbon) %>%
    data.frame() |>
    select(33) |>
    ggplot(aes(x = c)) + 
    geom_density(fill = 'gray') + 
    geom_vline(xintercept = 0) + 
    theme_bw() +
    theme(axis.text.y = element_blank())+
    labs(title = "Biomass effect", x = 'value'),
  as.matrix(fit_carbon) %>%
    data.frame() |>
    select(34) |>
    ggplot(aes(x = d)) + 
    geom_density(fill = 'gray') + 
    geom_vline(xintercept = 0) + 
    theme_bw() + 
    theme(axis.text.y = element_blank()) +
    labs(title = "Species richness effect", x = 'value'),
  heights = c(.7,.3)
)
#### >>> 3.2.2. Nitrogen model ####
#### >>>> 3.2.2.1 Predict original data ####
df = bind_cols(dat_stan_N.Ma[['p']], 
               dat_stan_N.Ma[['litterBM']],
               dat_stan_N.Ma[['litterSR']], 
               dat_stan_N.Ma[['decompC']])
df[,20] = NA
df = as.matrix(df)
for(i in 1: nrow(df)){
  df[i,20] = pred.decomp(newdata = df[i, 1:18], fit = fit_nitrogen)
}

#### >>>> 3.2.2.2. Bayesian R2 ####
R2.N = var(df[,20])/(var(df[,20]) + var(df[,20] - df[,19]))
R2.N

#### >>>> 3.2.2.3 Predicted vs measured ####
plot(df[,20]~df[,19])
abline(a = 0, b= 1)

#### >>>> 3.2.2.4 Residuals distribution ####
plot(sqrt(abs(df[,20]-df[,19]))~df[,19])

#### >>>> 3.2.2.5 Posterior distribution ####
ggpubr::ggarrange(
  as.matrix(fit_nitrogen) |>
    data.frame() |>
    select(1:16) |>
    pivot_longer(cols = 1:16) |>
    ggplot(aes(x = value, fill = name)) + 
    geom_density() + 
    geom_vline(xintercept = 0) + 
    facet_grid(rows = vars(name)) + 
    labs(title = 'Species-specific decomposition') + 
    theme_bw() +
    theme(axis.text.y = element_blank(), 
          legend.position = 'none',
          strip.text.y = element_blank()),
  
  as.matrix(fit_nitrogen) %>%
    data.frame() |>
    select(17:32) |>
    pivot_longer(cols = 1:16) %>%
    ggplot(aes(x = value, fill = name)) + 
    geom_density() + 
    geom_vline(xintercept = 0) + 
    facet_grid(rows = vars(name)) + 
    theme_bw()+
    labs(title = 'Interaction effect') + 
    theme(axis.text.y = element_blank(), 
          legend.position = 'none',
          strip.text.y = element_blank()),
  as.matrix(fit_nitrogen) %>%
    data.frame() |>
    select(33) |>
    ggplot(aes(x = c)) + 
    geom_density(fill = 'gray') + 
    geom_vline(xintercept = 0) + 
    theme_bw() +
    theme(axis.text.y = element_blank())+
    labs(title = "Biomass effect", x = 'value'),
  as.matrix(fit_nitrogen) %>%
    data.frame() |>
    select(34) |>
    ggplot(aes(x = d)) + 
    geom_density(fill = 'gray') + 
    geom_vline(xintercept = 0) + 
    theme_bw() + 
    theme(axis.text.y = element_blank()) +
    labs(title = "Species richness effect", x = 'value'),
  heights = c(.7,.3)
)

#### > 4. Saving for predictions later #### 
save.image(file = "1-2_fit-decomposition.RData")
#### END ####