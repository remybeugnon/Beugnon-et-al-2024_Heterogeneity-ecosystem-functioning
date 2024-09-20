# R version 4.3.1 (2023-06-16)
# Platform: aarch64-apple-darwin20 (64-bit)
# Running under: macOS Sonoma 14.6.1

# Packages
# tidyverse_2.0.0
# purrr_1.0.2

rm(list = ls())
set.seed(1603)

#### > 1. Functions ####
library(tidyverse)
library(purrr)

# Counting number of neighbors
nb.same.type = function(grid){
  # return the number of neighboors that have the same type
  sum.neighboors = function(c){
    i = c[2] # row
    j = c[1] # col
    x = c(i-1, i, i, i+1) 
    y = c(j, j+1, j-1, j)
    
    to.keep = x > 0 & x <= nrow & y > 0 & y <= ncol
    x = x[to.keep]
    y = y[to.keep]
    
    # move back to flat coordinates
    coords.flat = (x-1)*ncol + y
    return(sum(grid[(i-1)*ncol + j] == grid[coords.flat]))
  }
  mat.coords = expand.grid(1:ncol, 1:nrow)
  apply(mat.coords, MARGIN = 1, FUN = sum.neighboors)
}

# Measuring heterogeneity compared to hyperbolic law
expected.hypergeom = function(nb.n, nb.per.species){
  return(((nb.per.species-1)*nb.n)/(nrow*ncol-1))
}

# Heterogeneity measure
to.optimise = function(grid){
  nb.neighbours = c(c(2, rep(3, ncol-2), 2), rep(c(3,rep(4, ncol-2), 3), nrow-2), c(2, rep(3, ncol-2), 2)) 
  expectations.hypergeom = vapply(nb.neighbours, expected.hypergeom, FUN.VALUE = numeric(1), nb.per.species = nrow*ncol/length(unique(grid)))
  observed = nb.same.type(grid) 
  return(sum(observed-expectations.hypergeom))
}

# Random permutations of trees
permute = function(vect_forest, nb){
  for(i in 1:nb){
    permutation = sample(1:length(vect_forest), 2)
    a = vect_forest[permutation[1]]
    b = vect_forest[permutation[2]]
    vect_forest = 
      vect_forest %>% 
      replace(., permutation[1], b) %>%
      replace(., permutation[2], a)
  }
  return(vect_forest)
}

#### > 2. Setting up forest spatial structure #### 
ncol = 16
nrow = 16
grid = rep(LETTERS[1:8], ncol*nrow/8)

#### > 3. Design forests ####
#### >> 3.1. Fixed designs (block and rows)
designs.8 = {
  list(
    design.bloc = c(
      rep(c(rep('A',4),rep('B',4),rep('C',4),rep('D',4)),8),
      rep(c(rep('E',4),rep('F',4),rep('G',4),rep('H',4)),8)),
    
    design.minibloc = c(
      rep(c(rep('A',4),rep('B',4),rep('C',4),rep('D',4)),4),
      rep(c(rep('E',4),rep('F',4),rep('G',4),rep('H',4)),4),
      rep(c(rep('A',4),rep('B',4),rep('C',4),rep('D',4)),4),
      rep(c(rep('E',4),rep('F',4),rep('G',4),rep('H',4)),4)),
    
    design.2lines = c(
      rep("A", 32), rep("B", 32), 
      rep("C", 32), rep("D", 32), 
      rep("E", 32), rep("F", 32),
      rep("G", 32), rep("H", 32)),
    
    design.lines = c(
      rep("A", 16), rep("B", 16), 
      rep("C", 16), rep("D", 16), 
      rep("E", 16), rep("F", 16),
      rep("G", 16), rep("H", 16),
      rep("A", 16), rep("B", 16), 
      rep("C", 16), rep("D", 16), 
      rep("E", 16), rep("F", 16),
      rep("G", 16), rep("H", 16))
  )
}

designs.4 = {
  list(
    design.bloc = c(
      rep(c(rep('A',8),rep('B',8)),8),
      rep(c(rep('C',8),rep('D',8)),8)
    ),
    design.minibloc = c(
      rep(c(rep('A',4),rep('B',4),rep('C',4),rep('D',4)),4),
      rep(c(rep('C',4),rep('D',4),rep('A',4),rep('B',4)),4),
      rep(c(rep('A',4),rep('B',4),rep('C',4),rep('D',4)),4),
      rep(c(rep('C',4),rep('D',4),rep('A',4),rep('B',4)),4)
    ),
    design.2lines = 
      rep(c(rep('A',32),rep('B',32),rep('A',32),rep('B',32)), 2),
    
    design.lines = 
      rep(c(rep('A',16),rep('B',16),rep('A',16),rep('B',16)), 4)
  )
}

designs.2 = {
  list(
    design.bloc = 
      c(rep('A',16*8),rep('B',16*8)),
    
    design.minibloc = 
      c(
        rep(c(rep(c(rep('A',4),rep('B',4)),2)),4),
        rep(c(rep(c(rep('B',4),rep('A',4)),2)),4),
        rep(c(rep(c(rep('A',4),rep('B',4)),2)),4),
        rep(c(rep(c(rep('B',4),rep('A',4)),2)),4)
      ),
    
    design.2lines = 
      rep(c(rep('A',32),rep('B',32)), 4),
    
    design.lines = 
      rep(c(rep('A',16),rep('B',16)), 8)
  )
}

# Random permutation from block design
designs.8 = c(
  designs.8, 
  map(.x = 1:400, .f = function(x){
    forest = permute(designs.8[['design.bloc']], x)
    forest})
)

designs.4 = c(
  designs.4, 
  map(.x = 1:500, .f = function(x){
    forest = permute(designs.4[['design.bloc']], x)
    forest})
)

designs.2 = c(
  designs.2, 
  map(.x = 1:500, .f = function(x){
    forest = permute(designs.2[['design.bloc']], x)
    forest})
)

#### >> 3.2. Forest selections ####
forest.designs.8 = designs.8[c(1:4,15,25,45,50,75,110,140,305)]
forest.designs.4 = designs.4[c(1:4,12,25,47,67,105,155,500)]
forest.designs.2 = designs.2[c(1:4,12,26,47,75,114,175,400)]

#### >> 3.3. Heterogeneity calculation ####
heterogeneity.forest.8 =
  data.frame(
    forests.8 = 1:12,
    description = c(
      'Blocks', 'Miniblock', 'Two rows', 'One row',
      paste(c(15,25,44,56,75,110,155,305), 'permutations')
    ) %>% factor(., levels = c(
      'Blocks', 'Miniblock', 'Two rows', 'One row', 
      paste(c(15,25,44,56,75,110,155,305), 'permutations')
    )),
    homogeneity = map_dbl(.x = 1:12, 
                          .f = function(x){
                            to.optimise(forest.designs.8[[x]])})
  )

heterogeneity.forest.4 =
  data.frame(
    forests.4 = 1:11,
    description = c(
      'Blocks', 'Miniblock', 'Two rows', 'One row',
      paste(c(12,25,47,75,100,165,500), 'permutations')
    ) %>% factor(., levels = c(
      'Blocks', 'Miniblock', 'Two rows', 'One row',
      paste(c(12,25,47,75,100,165,500), 'permutations')
    )),
    homogeneity = map_dbl(.x = 1:11, 
                          .f = function(x){
                            to.optimise(forest.designs.4[[x]])})
  )

heterogeneity.forest.2 =
  data.frame(
    forests.2 = 1:9,
    description = c(
      'Blocks', 'Miniblock', 
      paste(c(12,26,47,75,114,165,400), 'permutations')
    ) %>% factor(., levels = c(
      'Blocks', 'Miniblock',
      paste(c(12,26,47,75,114,165,400), 'permutations')
    )),
    homogeneity = map_dbl(.x = 1:9, 
                          .f = function(x){
                            to.optimise(forest.designs.2[[x]])})
  )

# Plot heterogeneity
# plot(rep(1,9)~heterogeneity.forest.8$homogeneity[c(1,5:12)])
# plot(rep(1,8)~heterogeneity.forest.4$homogeneity[c(1,5:11)])
# plot(rep(1,6)~heterogeneity.forest.4$homogeneity[c(1,5:9)])

#### > 4. Saving ####
saveRDS(forest.designs.2, file = '2-1-1_2sp-forests.RDS')
saveRDS(forest.designs.4, file = '2-1-1_4sp-forests.RDS')
saveRDS(forest.designs.8, file = '2-1-1_8sp-forests.RDS')
#### END ####