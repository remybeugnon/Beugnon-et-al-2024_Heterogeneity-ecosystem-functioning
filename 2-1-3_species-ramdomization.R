# R version 4.3.1 (2023-06-16)
# Platform: aarch64-apple-darwin20 (64-bit)
# Running under: macOS Sonoma 14.6.1

# Packages
# tidyverse_2.0.0
rm(list = ls())
set.seed(1603)
library(tidyverse)

# Loading data
species.poll = read_rds(file = '2-1-2_forests-arrangements.RDS')

replace.species.1 = function(x, i){
  y = x |>
    str_replace_all(LETTERS[1], species.poll[[1]][i,1])
  return(y)
}

replace.species.2 = function(x, i){
  y = x |>
    str_replace_all(LETTERS[1], as.character(species.poll[[2]][i,1])) |>
    str_replace_all(LETTERS[2], as.character(species.poll[[2]][i,2])) 
  return(y)
}

replace.species.4 = function(x, i){
  y = x |>
    str_replace_all(LETTERS[1], as.character(species.poll[[3]][i,1])) |>
    str_replace_all(LETTERS[2], as.character(species.poll[[3]][i,2])) |>
    str_replace_all(LETTERS[3], as.character(species.poll[[3]][i,3])) |>
    str_replace_all(LETTERS[4], as.character(species.poll[[3]][i,4]))
  return(y)
}

replace.species.8 = function(x, i){
  y = x |>
    str_replace_all(LETTERS[1], as.character(species.poll[[4]][i,1])) |>
    str_replace_all(LETTERS[2], as.character(species.poll[[4]][i,2])) |>
    str_replace_all(LETTERS[3], as.character(species.poll[[4]][i,3])) |>
    str_replace_all(LETTERS[4], as.character(species.poll[[4]][i,4])) |>
    str_replace_all(LETTERS[5], as.character(species.poll[[4]][i,5])) |>
    str_replace_all(LETTERS[6], as.character(species.poll[[4]][i,6])) |>
    str_replace_all(LETTERS[7], as.character(species.poll[[4]][i,7])) |>
    str_replace_all(LETTERS[8], as.character(species.poll[[4]][i,8]))
  return(y)
}

# Adding the species permutations to the design 
forests.2 = list()
forests.2[['design']] = readRDS(file = '2-1-1_2sp-forests.RDS')
names(forests.2[["design"]]) = c('block', 'mini.block', 'two.lines', 'one.line',
                                 paste0(c(12,26,47,75,114,165,400), 'perms'))
forests.2[['species']] = vector(mode = 'list', length = length(forests.2[['design']]))
names(forests.2[["species"]]) = names(forests.2[["design"]])

for(i in 1:length(forests.2[['design']])){
  forests.2[["species"]][[i]] = vector(mode = 'list', length = nrow(species.poll[[2]]))
  for(j in 1:nrow(species.poll[[2]]))
    forests.2[['species']][[i]][[j]] = replace.species.2(forests.2[['design']][[i]],j)
}

forests.4 = list()
forests.4[['design']] = readRDS(file = '2-1-1_4sp-forests.RDS')
names(forests.4[["design"]]) = c('block', 'mini.block', 'two.lines', 'one.line',
                                 paste0(c(12,25,47,75,100,165,500), 'perms'))
forests.4[['species']] = vector(mode = 'list', length = length(forests.4[['design']]))
names(forests.4[["species"]]) = names(forests.4[["design"]])

for(i in 1:length(forests.4[['design']])){
  forests.4[["species"]][[i]] = vector(mode = 'list', length = nrow(species.poll[[3]]))
  for(j in 1:nrow(species.poll[[3]]))
    forests.4[['species']][[i]][[j]] = replace.species.4(forests.4[['design']][[i]],j)
}

forests.8 = list()
forests.8[['design']] = readRDS(file = '2-1-1_8sp-forests.RDS')
names(forests.8[["design"]]) = c('block', 'mini.block', 'two.lines', 'one.line',
                               paste0(c(15,25,44,56,75,110,155,305), 'perms'))
forests.8[['species']] = vector(mode = 'list', length = length(forests.8[['design']]))
names(forests.8[["species"]]) = names(forests.8[["design"]])

for(i in 1:length(forests.8[['design']])){
  forests.8[["species"]][[i]] = vector(mode = 'list', length = nrow(species.poll[[4]]))
  for(j in 1:nrow(species.poll[[4]]))
  forests.8[['species']][[i]][[j]] = replace.species.8(forests.8[['design']][[i]],j)
}

# Saving th data
saveRDS(forests.2, file = '2-1-3_species-randomization-2sp.RDS')
saveRDS(forests.4, file = '2-1-3_species-randomization-4sp.RDS')
saveRDS(forests.8, file = '2-1-3_species-randomization-8sp.RDS')
#### END ####