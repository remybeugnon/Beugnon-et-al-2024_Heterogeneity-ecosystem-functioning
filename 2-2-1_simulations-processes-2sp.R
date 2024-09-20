# R version 4.3.1 (2023-06-16)
# Platform: aarch64-apple-darwin20 (64-bit)
# Running under: macOS Sonoma 14.6.1

rm(list = ls())

# Setting start and stop for the server parallel calculation
start = 1
stop  = 2

set.seed(1603)

# Packages
# tidyverse_2.0.0

library(tidyverse)
source('1-3_prediction-function.R')

# Loading forest structure
forests = readRDS(file = '2-1-3_species-randomization-2sp.RDS')

# Simulations
for(i in start:stop){
    result.sim = data.frame(matrix(data = NA,
                                   nrow = length(forests$design),
                                   ncol = 14))
    colnames(result.sim) = c("design", "tree.biomass", "tree.biomass.sd",
                             "total.litter", "mean.litter",
                             "sd.litter", "mean.sp.rich.litter",
                             "sd.sp.rich.litter",   "mean.decomp.rate.C",
                             "sd.decomp.rate.C",    "mean.decomp.rate.N",
                             "sd.decomp.rate.N",    "total.decomp.C",
                             "total.decomp.N")
    
    for(j in 1:length(forests[['design']])){
      s.p = simul(forest = forests$species[[j]][[i]])
      # Saving raw data
      write.csv(s.p[['final_forest']],
                paste0('simul/2/simul-2sp_design-', j, '_rep-',i,'.csv'),
                row.names = F)
      result.sim[j,] = c(j, s.p[['summary']])
    }
    # Saving the summary
    write.csv(result.sim,
              paste0('simul/output-2sp_rep-',i,'.csv'),
              row.names = F)
}
#### END ####