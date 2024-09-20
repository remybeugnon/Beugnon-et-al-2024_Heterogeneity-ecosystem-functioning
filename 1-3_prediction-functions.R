# R version 4.3.1 (2023-06-16)
# Platform: aarch64-apple-darwin20 (64-bit)
# Running under: macOS Sonoma 14.6.1

#### > 1. Data handling ####
pars <- read.csv('df-tree-interaction-parameters.csv', 
                 nrows = 74, row.names = 1)
alpha_matrix <- matrix(pars[grep("^alpha", row.names(pars)), "mean"], 
                       nrow = 8, byrow = TRUE)
vec_beta_sp <- pars[grep("^beta", row.names(pars)), "mean"]
theta <- pars["theta", "mean"]
b <- pars["b", "mean"]

# - for reference, species names and their IDs (as we used them)
# ID  - species
# S1  - Castanea henryi
# S2  - Castanopsis sclerophylla
# S3  - Choerospondias axillaris
# S4  - Liquidambar formosana
# S5  - Nyssa sinensis
# S6  - Quercus serrata
# S7  - Sapindus saponaria
# S8  - Triadica sebifera

#### > 2. define functions ####
## identify direct neighbours
# - dat is a dataframe that has treeIDs, x, and y coordinates
# - ind is identifier for individual trees in plot
# - indID can be used if values in ind don't use the format "ind1"
# - wrap specifies if tree interactions should be considered with periodic boundary conditions (TRUE) or not (FALSE)
neigh.wrap_function <- function(dat, ind, indID = "ind", wrap = TRUE){
  ## define which are edge trees and need wrap around (if we wrap!)
  if(wrap){
    
    edge <- dat[dat$x %in% c(min(dat$x), max(dat$x)) | 
                  dat$y %in% c(min(dat$y), max(dat$y)), ]
    edge$x[edge$x == min(dat$x)] <- max(dat$x) + 1
    edge$x[edge$x == max(dat$x)] <- min(dat$x) - 1
    edge$y[edge$y == min(dat$y)] <- max(dat$y) + 1
    edge$y[edge$y == max(dat$y)] <- min(dat$y) - 1
    
    edge_xfix <- edge_yfix <- edge[paste(edge$x, edge$y) %in% paste(rep(range(edge$x), 2), 
                                                                    rep(range(edge$y), each = 2)), ]
    edge_xfix$x <- abs(edge_xfix$x - max(dat$x))
    edge_yfix$y <- abs(edge_yfix$y - max(dat$y))
    
    
    edge <- rbind(edge, edge_xfix, edge_yfix)
    
    
    tmp.dat <- rbind(cbind(dat, core = rep("core", nrow(dat))), cbind(edge, core = rep("edge", nrow(edge))))
    
  } else {
    
    tmp.dat <- cbind(dat, core = rep("core", nrow(dat)))
    
  }
  
  ## specify neighbouring tree IDs !
  tmp <- tmp.dat %>%
    {.[.[, indID] == ind & .$core == "core", ]} %>%
    {tmp.dat[tmp.dat$x %in% (.$x-1):(.$x+1) & 
               tmp.dat$y %in% (.$y-1):(.$y+1) &
               paste(tmp.dat$x, tmp.dat$y) != paste(.$x, .$y), ]} %>%
    {cbind(.[, indID], 1)} %>%
    as.data.frame()
  colnames(tmp) <- c(indID, "dist")
  
  
  ## prepare output
  out <- suppressMessages(left_join(dat[indID], tmp)[, 2])
  return(as.numeric(out))
}


## simulating function modified from Yu et al. 2023 -- DOI: 10.1111/ele.14338  
# - dat is the imput data.frame with columns 'x', 'y', and 'species' (s. above)
# - years define how many years of growth should be simulated; default is 7 (i.e. 2009-2016), which is the time frame for which we got our parameter estimates
# - start.biomass defines the initial value for biomass; we set it to 100, assuming equal starting conditions
# - output defines if you get an output for the final years ('final') or all years ('all'); the latter adds column year to the output 

simgrowth <- function(forest, years = 10, start.biomass = 100, output = "final"){
  dat <- data.frame(
    expand_grid(x = 1:16, 
                y = 1:16)) |> 
    mutate(species = as.numeric(forest))
  
  ## number of individuals:
  nind <- nrow(dat)
  
  ## unique IDs for each tree:
  dat$ind <- paste0("ind", 1:nrow(dat))
  
  ## neighbor identifier
  df_neigh <- sapply(dat$ind, 
                     function(x) neigh.wrap_function(dat, x))
  colnames(df_neigh) <- dat$ind
  rownames(df_neigh) <- dat$ind
  
  ## beta vector (ind)
  vec_beta_ind <- vec_beta_sp[dat$species]
  
  ## alpha matrix (ind x ind)
  df_alpha_ind <- matrix(ncol = nind, nrow = nind)
  for(k in 1:ncol(df_alpha_ind)){
    df_alpha_ind[, k] <- alpha_matrix[, dat$species[k]][dat$species]
  }
  colnames(df_alpha_ind) <- dat$ind
  rownames(df_alpha_ind) <- dat$ind
  
  ## starting densities:
  biomass_df <- as.data.frame(matrix(rep(start.biomass, nind), nrow = 1))
  colnames(biomass_df) <- dat$ind
  
  ## simulate growth over years 
  for(k in 1:(years-1)){
    
    biomass <- as.numeric(biomass_df[k, ])
    
    ## calculate intrinsic growth term
    intrinsic_growth <- vec_beta_ind * biomass^theta
    
    
    ## calculate interaction term
    interaction_sum <- (df_alpha_ind * df_neigh) %>% apply(1, function(x){x * biomass^b}) %>% colSums(na.rm = TRUE)
    
    
    ## calculate new biomass
    biomass_df[k+1, ] <- biomass + intrinsic_growth + interaction_sum
    
    
    ## set negative biomasses to zero (rare!) & and assure dead trees stay dead
    biomass_df[k+1, biomass_df[k+1, ] <= 0] <- 0.0
    biomass_df[k+1, biomass_df[k, ] <= 0] <- 0.0
    
  }
  
  
  ## create output
  biomass_df$year <- 1:years
  
  out <- biomass_df %>% 
    pivot_longer(cols = starts_with("ind"), names_to = "ID_focal", values_to = "biomass") %>%
    as.data.frame()
  
  out$x <- rep(dat$x, years)
  out$y <- rep(dat$y, years)
  out$species <- rep(dat$species, years)
  
  
  if(output == "final"){
    return(out[out$year == years, c("x", "y", "species", "biomass")])
  } else if(output == "all"){
    return(out[c("x", "y", "species", "biomass", "year")])
  } else {
    warning("please specify output as 'all' or 'final' for either getting all years or only the last one")
  }
  
}

#### > 3. Litterfall and decomposition prediction ####
simul = function(forest, restrict.dispersal = T) {
  # Define grid
  point.grid = {
    expand.grid(x = seq(1, 16, 0.1),
                y = seq(1, 16, 0.1)) |>
      # add one column per species for distance, biomass and litter
      mutate(bio.S1 = 0) |>
      mutate(bio.S2 = 0) |>
      mutate(bio.S3 = 0) |>
      mutate(bio.S4 = 0) |>
      mutate(bio.S5 = 0) |>
      mutate(bio.S6 = 0) |>
      mutate(bio.S7 = 0) |>
      mutate(bio.S8 = 0) |>
      mutate(bio.S9 = 0) |>
      mutate(bio.S10 = 0) |>
      mutate(bio.S11 = 0) |>
      mutate(bio.S12 = 0) |>
      mutate(dist.S1 = 0) |>
      mutate(dist.S2 = 0) |>
      mutate(dist.S3 = 0) |>
      mutate(dist.S4 = 0) |>
      mutate(dist.S5 = 0) |>
      mutate(dist.S6 = 0) |>
      mutate(dist.S7 = 0) |>
      mutate(dist.S8 = 0) |>
      mutate(dist.S9 = 0) |>
      mutate(dist.S10 = 0) |>
      mutate(dist.S11 = 0) |>
      mutate(dist.S12 = 0) |>
      mutate(biodist.S1 = 0) |>
      mutate(biodist.S2 = 0) |>
      mutate(biodist.S3 = 0) |>
      mutate(biodist.S4 = 0) |>
      mutate(biodist.S5 = 0) |>
      mutate(biodist.S6 = 0) |>
      mutate(biodist.S7 = 0) |>
      mutate(biodist.S8 = 0) |>
      mutate(biodist.S9 = 0) |>
      mutate(biodist.S10 = 0) |>
      mutate(biodist.S11 = 0) |>
      mutate(biodist.S12 = 0) |>
      mutate(litter.S1 = 0) |>
      mutate(litter.S2 = 0) |>
      mutate(litter.S3 = 0) |>
      mutate(litter.S4 = 0) |>
      mutate(litter.S5 = 0) |>
      mutate(litter.S6 = 0) |>
      mutate(litter.S7 = 0) |>
      mutate(litter.S8 = 0) |>
      mutate(litter.S9 = 0) |>
      mutate(litter.S10 = 0) |>
      mutate(litter.S11 = 0) |>
      mutate(litter.S12 = 0)
  }
  
  # Calculate distance from trees
  forest.biomass = simgrowth(forest, years = 10) |>
    mutate(biomass = biomass/1000000)
  
  for (i in 1:nrow(point.grid)) {
    for (j in 1:nrow(forest.biomass)) {
      dist = sqrt((point.grid$x[i] - forest.biomass$x[j]) ^ 2 +
                    (point.grid$y[i] - forest.biomass$y[j]) ^ 2)
      if (restrict.dispersal == T) {
        if (dist < sqrt(1.5 ^ 2 + 1) & dist > 0) {
          sp = forest.biomass$species[j]
          point.grid[i, paste0('bio.S', sp)] = point.grid[i, paste0('bio.S', sp)] + forest.biomass$biomass[j] 
          point.grid[i, paste0('dist.S', sp)] = point.grid[i, paste0('dist.S', sp)] + (1 / dist)
          point.grid[i, paste0('biodist.S', sp)] = point.grid[i, paste0('biodist.S', sp)] +  (forest.biomass$biomass[j] / dist)
        }
      }else{
        if (dist > 0) {
          sp = forest.biomass$species[j]
          point.grid[i, paste0('bio.S', sp)] = point.grid[i, paste0('bio.S', sp)] + forest.biomass$biomass[j] 
          point.grid[i, paste0('dist.S', sp)] = point.grid[i, paste0('dist.S', sp)] + (1 / dist)
          point.grid[i, paste0('biodist.S', sp)] = point.grid[i, paste0('biodist.S', sp)] +  (forest.biomass$biomass[j] / dist)
        }
      }
    }
  }
  
  #### > 4. litter distribution ####
  #### >> 4.1. Extraction model fit ####
  load(file = "fit-litter.RData")
  
  parameters.b =
    as.matrix(fit) |>
    as.data.frame() |>
    select(paste0('bio[', 1:12, ']')) |>
    colMeans()
  
  parameters.d =
    as.matrix(fit) |>
    as.data.frame() |>
    select(paste0('d[', 1:12, ']')) |>
    colMeans()
  
  parameters.db =
    as.matrix(fit) |>
    as.data.frame() |>
    select(paste0('db[', 1:12, ']')) |>
    colMeans()
  
  #### >> 4.2. Estimate species specific fall ####
  for (i in 1:nrow(point.grid)) {
    for (sp in 1:12) {
      est =
        parameters.b[sp] * point.grid[i, paste0('bio.S', sp)] +
        parameters.d[sp] * point.grid[i, paste0('dist.S', sp)] +
        parameters.db[sp] * point.grid[i, paste0('biodist.S', sp)]
      if (est > 0) {
        point.grid[i, paste0('litter.S', sp)] = est
      }
    }
  }
  
  #### >> 4.3. Aggregate indices ####
  point.grid$litter.total = rowSums(point.grid[, grep('litter.S', colnames(point.grid))])
  point.grid$litter.sp.rich = rowSums(point.grid[, grep('litter.S', colnames(point.grid))] !=
                                        0)
  
  #### > 5. Predicting decomposition ####
  #### >> 5.1. Extraction model fit ####
  load(file = "fit-decomposition.RData")
  
  parameters.decom.C = fit_carbon |>
    as.matrix() |>
    colMeans()
  
  parameters.decom.N = fit_nitrogen |>
    as.matrix() |>
    colMeans()
  
  #### >> 5.2. Predict C and N loss ####
  point.grid =
    bind_cols(point.grid,
              point.grid[,
                         grep('litter.S',
                              colnames(point.grid))] /
                point.grid$litter.total)
  
  point.grid$decomp.C = NA
  point.grid$decomp.N = NA
  
  for (i in 1:nrow(point.grid)) {
    new.data = point.grid |>
      as.matrix()
    newdata = new.data[i, c(53:64, 51, 52)]
    par.C = fit_carbon |>
      as.matrix() |>
      colMeans()
    par.N = fit_nitrogen |>
      as.matrix() |>
      colMeans()
    
    pred.c =
      sum(par.C[1:12] * newdata[1:12]) + # identity effect
      sum(par.C[17:28] * newdata[1:12] * (1 - newdata[1:12])) + # interaction effect
      par.C[33] * newdata[13] + # biomass effect
      par.C[34] * newdata[14]
    
    if (pred.c > 100) {
      point.grid$decomp.C[i] =  100
    } else if (pred.c < 0) {
      point.grid$decomp.C[i] =  0
    } else{
      point.grid$decomp.C[i] = pred.c
    }
    
    pred.n =
      sum(par.N[1:12] * newdata[1:12]) + # identity effect
      sum(par.N[17:28] * newdata[1:12] * (1 - newdata[1:12])) + # interaction effect
      par.N[33] * newdata[13] + # biomass effect
      par.N[34] * newdata[14]
    
    if (pred.n > 100) {
      point.grid$decomp.N[i] =  100
    } else if (pred.n < 0) {
      point.grid$decomp.N[i] =  0
    } else{
      point.grid$decomp.N[i] = pred.n
    }
  }
  #### > 6. Export ####
  return(list(
    final_forest = point.grid,
    summary = c(tree.biomass = sum(forest.biomass$biomass),
                tree.biomass.sd = sd(forest.biomass$biomass),
                unlist(summarise(
                  point.grid,
                  total.litter = sum(litter.total),
                  mean.litter = mean(litter.total),
                  sd.litter = sd(litter.total),
                  mean.sp.rich.litter = mean(litter.sp.rich),
                  sd.sp.rich.litter = sd(litter.sp.rich),
                  mean.decomp.rate.C = mean(decomp.C),
                  sd.decomp.rate.C = sd(decomp.C),
                  mean.decomp.rate.N = mean(decomp.N),
                  sd.decomp.rate.N = sd(decomp.N),
                  total.decomp.C = sum(decomp.C * litter.total * 0.5  / 100),
                  total.decomp.N = sum(decomp.N * litter.total * 0.04 / 100)
                )))
  )
  )
}
#### END ####