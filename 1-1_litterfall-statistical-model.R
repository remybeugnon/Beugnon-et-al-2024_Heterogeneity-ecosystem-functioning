# R version 4.3.1 (2023-06-16)
# Platform: aarch64-apple-darwin20 (64-bit)
# Running under: macOS Sonoma 14.6.1

# Packages
# tidyverse_2.0.0
# rstan_2.33.1.9000 
# coda_0.19-4
# loo_2.6.0

rm(list = ls())

#### > 1. Initialization ####
#### >> 1.1. Packages ####
libs <- c(
  'tidyverse',
  'rstan', 'bayesplot', 
  'coda', 'loo' 
)
invisible(lapply(libs, library, character.only = T))

#### >> 1.2. Datasets ####
df = read.csv('df-litterfall.csv') %>% 
  select(-dist, -biomass, -n.rich)

df.tree = read.csv('df-tree.csv') %>%
  select(TSP, species, biomass, dist, db)

df = df %>%
  left_join(., 
            df.tree, 
            by = c('TSP','species'))

# Coding species names
df$species = df$species %>%
  str_replace_all( . , pattern = "Castanea henryi", 'A') %>%
  str_replace_all( . , pattern = "Castanopsis sclerophylla", 'B') %>%
  str_replace_all( . , pattern = "Choerospondias axillaris", 'C') %>%
  str_replace_all( . , pattern = "Liquidambar formosana", 'D') %>%
  str_replace_all( . , pattern = "Nyssa sinensis", 'E') %>%
  str_replace_all( . , pattern = "Quercus serrata", 'F') %>%
  str_replace_all( . , pattern = "Sapindus mukorossi", 'G') %>%
  str_replace_all( . , pattern = "Sapium sebiferum", 'H') %>%
  str_replace_all( . , pattern = "Quercus fabri", 'I') %>%
  str_replace_all( . , pattern = "Lithocarpus glaber", 'J') %>%
  str_replace_all( . , pattern = "Koelreuteria bipinnata", 'K') %>%
  str_replace_all( . , pattern = "Cyclobalanopsis glauca", 'L') 

#### >> 1.3. Data handling ####
df.1 =
  df %>%
  select(TSP, 
         litter = litter.biomass.area, 
         species) %>%
  filter(!is.na(species)) %>%
  group_by(TSP, species) %>%
  summarise(litter = sum(litter)) %>%
  pivot_wider(.,
              id_cols = TSP,
              names_from = species, 
              values_from = litter, 
              values_fill = 0,
              names_prefix = "litter_",
              names_sort = T) %>%
  left_join(.,
            df %>%
              select(TSP, 
                     biomass = biomass, 
                     species) %>%
              filter(!is.na(species)) %>%
              group_by(TSP, species) %>%
              summarise(biomass = sum(biomass)) %>%
              pivot_wider(.,
                          id_cols = TSP,
                          names_from = species, 
                          values_from = biomass, 
                          values_fill = 0,
                          names_prefix = "biomass_",
                          names_sort = T)
  ) %>%
  left_join(.,
            df %>%
              select(TSP, 
                     distance = dist, 
                     species) %>%
              filter(!is.na(species)) %>%
              group_by(TSP, species) %>%
              summarise(distance = sum(distance)) %>%
              pivot_wider(.,
                          id_cols = TSP,
                          names_from = species, 
                          values_from = distance, 
                          values_fill = 0,
                          names_prefix = "dist_",
                          names_sort = T)
  ) %>%
  left_join(.,
            df %>%
              select(TSP, 
                     biodist = db, 
                     species) %>%
              filter(!is.na(species)) %>%
              group_by(TSP, species) %>%
              summarise(biodist = sum(biodist)) %>%
              pivot_wider(.,
                          id_cols = TSP,
                          names_from = species, 
                          values_from = biodist, 
                          values_fill = 0,
                          names_prefix = "biodist_",
                          names_sort = T)
  )

#### > 2. Stan model ####
#### >> 2.1 Model ####
stan_code = '
data{
  int n; // observations per species
  int m; // species number
  matrix[m,n] litter; // Species specific litter
  matrix[m,n] biomass;  // Species specific biomass
  matrix[m,n] dist;  // Species specific distance
  matrix[m,n] biodist;  // Species specific distance
}

parameters{
  vector<lower=0>[m] bio; // Species specific biomass parameter
  vector<lower=0>[m] d; // Species specific distance parameter
  vector<lower=0>[m] db;  // Species specific interraction parameter
  vector<lower=0>[m] sigma;  // Species specific sigma
}

model{
  // Priors
  matrix[m,n] mu;  // Species specific mu 
  bio ~ normal(0,10); // Species specific biomass parameter prior
  d ~ normal(0,10);  // Species specific distance parameter
  db ~ normal(0,10);  // Species specific interaction parameter
  sigma ~ normal(0,10); // Species specific sigma

  // Likelihood
  for(i in 1:n){
    for(j in 1:m){
      mu[j,i] = bio[j] * biomass[j,i] + d[j] * dist[j,i] + db[j] * biodist[j,i] ; //Species specific litterfall
      litter[j,i] ~ normal(mu[j,i], sigma[j]); 
    }
  }
}
'

stan_model = stan_model(model_code = stan_code)

data = list(n = nrow(df.1),
            m = 12,
            litter = rbind(df.1$litter_A,
                           df.1$litter_B,
                           df.1$litter_C,
                           df.1$litter_D,
                           df.1$litter_E,
                           df.1$litter_F,
                           df.1$litter_G,
                           df.1$litter_H,
                           df.1$litter_I,
                           df.1$litter_J,
                           df.1$litter_K,
                           df.1$litter_L
            ),
            biomass = rbind(df.1$biomass_A,
                            df.1$biomass_B,
                            df.1$biomass_C,
                            df.1$biomass_D,
                            df.1$biomass_E,
                            df.1$biomass_F,
                            df.1$biomass_G,
                            df.1$biomass_H,
                            df.1$biomass_I,
                            df.1$biomass_J,
                            df.1$biomass_K,
                            df.1$biomass_L
            ),
            dist = rbind(df.1$dist_A,
                         df.1$dist_B,
                         df.1$dist_C,
                         df.1$dist_D,
                         df.1$dist_E,
                         df.1$dist_F,
                         df.1$dist_G,
                         df.1$dist_H,
                         df.1$dist_I,
                         df.1$dist_J,
                         df.1$dist_K,
                         df.1$dist_L
            ),
            biodist = rbind(df.1$biodist_A,
                            df.1$biodist_B,
                            df.1$biodist_C,
                            df.1$biodist_D,
                            df.1$biodist_E,
                            df.1$biodist_F,
                            df.1$biodist_G,
                            df.1$biodist_H,
                            df.1$biodist_I,
                            df.1$biodist_J,
                            df.1$biodist_K,
                            df.1$biodist_L
            )
)

#### >> 2.2. Fit model ####
fit = sampling(stan_model, 
               data = data,
               iter = 3000,
               warmup = 1000,
               chains = 4,
               cores = 4)

#### >> 2.3. Output ####
print(fit)

check_hmc_diagnostics(fit)

posterior = As.mcmc.list(fit)
plot(posterior)


# Extract data posterior
# Biomass
post.bio = as.matrix(fit) %>%
  as.data.frame() %>%
  select(c(paste0('bio[', 1:12,']'))) 

ggplot(data = post.bio |> pivot_longer(1:12), 
       aes(x = value, color = name)) + 
  geom_density()

post.bio %>% colMeans()

# Distance
post.d = as.matrix(fit) %>%
  as.data.frame() %>%
  select(c(paste0('d[', 1:12,']'))) 

ggplot(data = post.d |> pivot_longer(1:12), 
       aes(x = value, color = name)) + 
  geom_density()

post.d %>% colMeans()

# Biomass/Distance
post.db = as.matrix(fit) %>%
  as.data.frame() %>%
  select(c(paste0('db[', 1:12,']'))) 

post.db %>% colMeans()

ggplot(data = post.db |> pivot_longer(1:12), 
       aes(x = value, color = name)) + 
  geom_density()

post.sigma = as.matrix(fit) %>%
  as.data.frame() %>%
  select(c(paste0('sigma[', 1:12,']'))) 

#### >> 2.4. Model checking ####
n.post = nrow(post.bio)
mu = matrix(NA, nrow=data$m, data$n) # auxiliary variable for predictions
R2 = matrix(NA, nrow=n.post, ncol=data$m) # species-level R2
R2.all = rep(NA, n.post) # total R2
LL = array(NA, dim=c(n.post, data$m, data$n)) # log-likelihood 

data$litter.all = as.vector(data$litter)

for(k in 1:n.post){
  for(j in 1:data$m){
    # predictions vectorised over all n observations i
    mu[j, ] =
      post.bio[k,j] * data$biomass[j, ] + 
      post.d[k,j] * data$dist[j, ] + 
      post.db[k,j] * data$biomass[j, ] * data$dist[j, ]  
    # log-likelihood values
    LL[k,j, ] = dnorm(data$litter[j,], mean=mu[j, ], sd=post.sigma[k,j], log=TRUE) 
    # Bayesian R2 # https://doi.org/10.1080/00031305.2018.1549100
    R2[k,j] = var(mu[j, ]) / (var(mu[j, ])+var(mu[j, ]-data$litter[j, ])) 
  }
  # Bayesian R2 for joint dataset
  mu.all = as.vector(mu)
  R2.all[k] = var(mu.all) / (var(mu.all)+var(mu.all-data$litter.all)) 
}

# R2 estimates
colMeans(R2) # species-level
mean(R2.all) # total
# LOO-CV

ll = matrix(NA, nrow=n.post, ncol=data$m*data$n) # convert LL-array to "long" format
for(k in 1:n.post){
  ll[k, ] = as.vector(LL[k, , ])
}

loo.1 = loo(ll, cores=4)
loo.1

save.image(file = "1-1_fit-litterfall.RData")

# Check prediction
d.obs = 
  df.1 %>%
  select(TSP,
         litter_A, litter_B, litter_C,
         litter_D, litter_E, litter_F,
         litter_G, litter_H, litter_I,
         litter_J, litter_K, litter_L) %>%
  pivot_longer(cols = 2:13)

d.pred = 
  mu %>%
  t() %>%
  data.frame() %>%
  mutate(TSP = df.1$TSP)

colnames(d.pred) = c(paste0('litter_', LETTERS[1:12]),'TSP')

d.predi = d.pred %>%
  pivot_longer(cols = 1:12)

df.pred = 
  left_join(d.obs, 
            d.predi, 
            by = c('TSP','name')) |> 
  left_join(df, by = 'TSP')

ggplot(data = df.pred, 
       aes(x = value.y, y = value.x, color = name)) + 
  geom_point() +
  geom_abline(intercept = 0, slope = 1) +
  geom_abline(intercept = 0, slope = 0.80, lty = 2) +
  geom_abline(intercept = 0, slope = 1.20, lty = 2) +
  geom_abline(intercept = 0, slope = 0.50, lty = 3) +
  geom_abline(intercept = 0, slope = 1.50, lty = 3) +
  facet_wrap(vars(name)) +
  theme_bw() +
  coord_fixed()

ggplot(data = df.pred, 
       aes(x = value.y, y = value.x, color = name)) + 
  geom_point() +
  geom_abline(intercept = 0, slope = 0.80) +
  geom_abline(intercept = 0, slope = 1) +
  geom_abline(intercept = 0, slope = 1.20) +
  theme_bw() +
  coord_fixed()

ggplot(data = df.pred, 
       aes(x = spe.rich, y = value.x)) + 
  geom_jitter(aes(x = spe.rich, y = value.x)) +
  geom_smooth(method = 'loess', se = F, aes(x = spe.rich, y = value.x)) +
  geom_jitter(aes(x = spe.rich, y = value.y, color = name), alpha = .5, color = 'red') +
  geom_smooth(method = 'loess', se = F, aes(x = spe.rich, y = value.y), color = 'red', alpha = .5) + 
  facet_wrap(vars(name),scales = 'free') +
  theme_bw()

ggplot(data = df, aes(x = log(spe.rich), y = spe.rich * biomass/(12))) + 
  geom_jitter() + 
  geom_smooth(method = 'lm') + 
  facet_wrap(vars(species), scales = 'free')
#### END ####