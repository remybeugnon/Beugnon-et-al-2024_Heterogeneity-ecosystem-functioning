# R version 4.3.1 (2023-06-16)
# Platform: aarch64-apple-darwin20 (64-bit)
# Running under: macOS Sonoma 14.6.1
rm(list = ls())
set.seed(1603)

#### > 1. Packages and functions ####
library(tidyverse) # tidyverse_2.0.0
library(purrr) # purrr_1.0.2 

# Function
ncol = 16
nrow = 16
nb_species = 8
grid = rep(LETTERS[1:8], ncol*nrow/8)
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
expected.hypergeom = function(nb.n, nb.per.species){
  return(((nb.per.species-1)*nb.n)/(nrow*ncol-1))
}
to.optimise = function(grid, expected){
  nb.neighbours = c(c(2, rep(3, ncol-2), 2), rep(c(3,rep(4, ncol-2), 3), nrow-2), c(2, rep(3, ncol-2), 2)) 
  expectations.hypergeom = vapply(nb.neighbours, expected.hypergeom, FUN.VALUE = numeric(1), nb.per.species = nrow*ncol/9) # length(unique(grid)))
  observed = nb.same.type(grid) 
  return(sum(observed-expectations.hypergeom))
}

#### > 2. Loading data ####
# Forrest structure
forests.2 = readRDS(file = '2-8-simulations-V2-8sp-based/2-1-3_species-randomization-2sp.RDS')
forests.4 = readRDS(file = '2-8-simulations-V2-8sp-based/2-1-3_species-randomization-4sp.RDS')
forests.8 = readRDS(file = '2-8-simulations-V2-8sp-based/2-1-3_species-randomization-8sp.RDS')

# 2 species mixtures
list.records.2 = 
  list.files('2-8-simulations-V2-8sp-based/simul/2/')
list.records.2 = list.records.2[grep('output', list.records.2)]
x = list.records.2[1]
df.2sp = map_df(.x = list.records.2[-1], 
            function(x){
              data.frame(
                rep = x, 
                read.csv(paste0('2-8-simulations-V2-8sp-based/simul/2/',x)))
            }) 
name.short = colnames(df.2sp)
colnames(df.2sp) = c(
  'Replicate','Design' , 'Total tree biomass (T)', 'Tree biomass SD (T)', 
  'Total litterfall (g)', 'Litterfall (g/dm2)', 'Litterfall SD (g/dm2)',
  'Litter spe. richness (#/dm2)','Litter species richness SD (#/dm2)',
  'C decomposition rate (%)', 'C decomposition rate SD (%)',
  'N decomposition rate (%)', 'N decomposition rate SD (%)',
  'C loss (g/dm2)','N loss (g/dm2)'
)

df.2sp$Replicate = 
  df.2sp$Replicate |> 
  str_remove_all('output-2sp_rep-') |>
  str_remove_all('output-4sp_rep-') |>
  str_remove_all('output-8sp_rep-') |> 
  str_remove_all('.csv') |> 
  as.numeric()

df.2sp$heterogeneity = NA
for(i in 1:length(forests.2$design)){
  df.2sp$heterogeneity[df.2sp$Design == i] = -(to.optimise(forests.2$design[[i]])+1)
}

# 4 species mixtures
list.records.4 = 
  list.files('2-8-simulations-V2-8sp-based/simul/4/')
list.records.4 = list.records.4[grep('output', list.records.4)]
x = list.records.4[1]
df.4sp = map_df(.x = list.records.4[-1], 
                function(x){
                  data.frame(
                    rep = x, 
                    read.csv(paste0('2-8-simulations-V2-8sp-based/simul/4/',x)))
                }) 
name.short = colnames(df.4sp)
colnames(df.4sp) = c(
  'Replicate','Design' , 'Total tree biomass (T)', 'Tree biomass SD (T)', 
  'Total litterfall (g)', 'Litterfall (g/dm2)', 'Litterfall SD (g/dm2)',
  'Litter spe. richness (#/dm2)','Litter species richness SD (#/dm2)',
  'C decomposition rate (%)', 'C decomposition rate SD (%)',
  'N decomposition rate (%)', 'N decomposition rate SD (%)',
  'C loss (g/dm2)','N loss (g/dm2)'
)

df.4sp$Replicate = 
  df.4sp$Replicate |> 
  str_remove_all('output-2sp_rep-') |>
  str_remove_all('output-4sp_rep-') |>
  str_remove_all('output-8sp_rep-') |> 
  str_remove_all('.csv') |> 
  as.numeric()
df.4sp$heterogeneity = NA
for(i in 1:length(forests.2$design)){
  df.4sp$heterogeneity[df.4sp$Design == i] = -(to.optimise(forests.4$design[[i]])+1)
}

# 8 species mixtures
list.records.8 = 
  list.files('2-8-simulations-V2-8sp-based/simul/8/')
list.records.8 = list.records.8[grep('output', list.records.8)]
x = list.records.8[1]
df.8sp = map_df(.x = list.records.8[-1], 
                function(x){
                  data.frame(
                    rep = x, 
                    read.csv(paste0('2-8-simulations-V2-8sp-based/simul/8/',x)))
                }) 
name.short = colnames(df.8sp)
colnames(df.8sp) = c(
  'Replicate','Design' , 'Total tree biomass (T)', 'Tree biomass SD (T)', 
  'Total litterfall (g)', 'Litterfall (g/dm2)', 'Litterfall SD (g/dm2)',
  'Litter spe. richness (#/dm2)','Litter species richness SD (#/dm2)',
  'C decomposition rate (%)', 'C decomposition rate SD (%)',
  'N decomposition rate (%)', 'N decomposition rate SD (%)',
  'C loss (g/dm2)','N loss (g/dm2)'
)

df.8sp$Replicate = 
  df.8sp$Replicate |> 
  str_remove_all('output-2sp_rep-') |>
  str_remove_all('output-4sp_rep-') |>
  str_remove_all('output-8sp_rep-') |> 
  str_remove_all('.csv') |> 
  as.numeric()

df.8sp$heterogeneity = NA
for(i in 1:length(forests.8$design)){
  df.8sp$heterogeneity[df.8sp$Design == i] = -(to.optimise(forests.8$design[[i]])+1)
}

# Combining results
df = 
  df.2sp |> 
  mutate(sp = 2) |> 
  add_row(df.4sp |> 
            mutate(sp = 4)) |> 
  add_row(df.8sp |> 
            mutate(sp = 8))

df = df |> filter(!(sp == 8 & Design == 11))
df$Design[df$Design == 12] = 11

df$`C loss (g/dm2)` =
  df$`C loss (g/dm2)`/(17*17*100)
df$`N loss (g/dm2)` =
  df$`N loss (g/dm2)`/(17*17*100)

#### > 3. Plot heterogeneity by diversity levels ####
df$Design = as.factor(df$Design)

#### >> 3.1. Biomass ####
p.biomass = 
  ggplot() + 
  geom_smooth(data = df |> 
                select(Design, sp, 
                       `Total tree biomass (T)`, `Tree biomass SD (T)`) |> 
                pivot_longer(cols = 3:4),
              aes(x = sp, 
                  y = value, 
                  color = Design), method = 'loess',
              alpha = .1, se = T) + 
  scale_x_continuous(trans = 'log2', breaks = c(2,4,8)) +
  scale_color_manual(values = c(
    '1' = 'darkred',
    '2' = 'red', 
    '3' = 'darkblue',
    '4' = 'blue', '5' = 'gray70',
    '6' = 'gray60','7' = 'gray50','8' = 'gray40',
    '9' = 'gray30','10' = 'gray20','11' = 'black',
    '12' = 'black')) +
  facet_wrap(vars(name), scales = 'free_y') + 
  labs(x = '', 
       y = "Value", 
       title = "A. Tree biomass") +
  theme_bw() + 
  theme(legend.position = 'none')
p.biomass

#### >> 3.2. Litterfall ####
p.litter = 
  ggplot() + 
  geom_smooth(data = df |> 
                select(Design, sp, 
                       `Litterfall (g/dm2)`, 
                       `Litterfall SD (g/dm2)`, 
                       `Litter spe. richness (#/dm2)`) |> 
                pivot_longer(cols = 3:5) |> 
                mutate(name = factor(name, 
                                     levels = c('Litterfall (g/dm2)', 
                                                'Litterfall SD (g/dm2)', 
                                                'Litter spe. richness (#/dm2)'))),
              aes(x = sp, 
                  y = value, 
                  color = Design), method = 'loess',
              alpha = .1, se = T) + 
  scale_x_continuous(trans = 'log2', breaks = c(2,4,8)) +
  scale_color_manual(values = c(
    '1' = 'darkred',
    '2' = 'red', 
    '3' = 'darkblue',
    '4' = 'blue', '5' = 'gray70',
    '6' = 'gray60','7' = 'gray50','8' = 'gray40',
    '9' = 'gray30','10' = 'gray20','11' = 'black',
    '12' = 'black')) +
  labs(x = '', 
       y = "Value", 
       title = "B. Litterfall") +
  facet_wrap(vars(name), scales = 'free_y') + 
  theme_bw() + 
  theme(legend.position = 'none')

p.litter

#### >> 3.3. Decomposition ####
p.decomp = 
  ggplot() + 
  geom_smooth(data = df |> 
                select(Design, sp, 
                       `C decomposition rate (%)`, 
                       `C decomposition rate SD (%)`, 
                       `C loss (g/dm2)`) |> 
                pivot_longer(cols = 3:5),
              aes(x = sp, 
                  y = value, 
                  color = Design), method = 'loess',
              alpha = .1, se = T) + 
  scale_x_continuous(trans = 'log2', breaks = c(2,4,8)) +
  scale_color_manual(values = c(
    '1' = 'darkred',
    '2' = 'red', 
    '3' = 'darkblue',
    '4' = 'blue', '5' = 'gray70',
    '6' = 'gray60','7' = 'gray50','8' = 'gray40',
    '9' = 'gray30','10' = 'gray20','11' = 'black',
    '12' = 'black')) +
  facet_wrap(vars(name), scales = 'free_y') + 
  labs(x = 'Species richness', 
       y = "Value", 
       title = "C. Decomposition") +
  theme_bw() + 
  theme(legend.position = 'none')

#### >> 3.4. Saving the figure ####
Fig.3 = 
  ggpubr::ggarrange(
    ggpubr::ggarrange(p.biomass, 
                      ggplot(data=NULL) + theme_void(), 
                      nrow = 1, widths = c(2,1)),
    p.litter,
    p.decomp, 
    nrow = 3,
    align = 'hv')

Fig.3

ggsave(plot = Fig.3, 
       filename = '2-8-simulations-V2-8sp-based/3-2_Fig3.png',
       height = 19, width = 20, units = 'cm')

ggsave(plot = Fig.3, 
       filename = '2-8-simulations-V2-8sp-based/3-2_Fig3.svg',
       height = 20, width = 20, units = 'cm')

#### > 4. Statistical models ####
df = df |> filter(!(Design %in% c(2:4)))
df$sp = log2(df$sp)

mod = lm(data = df, 
   formula = 'log(`Total tree biomass (T)`)~sp*heterogeneity')
mod |> summary()

mod = lm(data = df, 
         formula = '`Tree biomass SD (T)`~sp*heterogeneity')
mod |> summary()

mod = lm(data = df, 
         formula = 'log(`Litterfall (g/dm2)`)~sp*heterogeneity')
mod |> summary()

mod = lm(data = df, 
         formula = '`Litterfall SD (g/dm2)`~sp*heterogeneity')
mod |> summary()

mod = lm(data = df, 
         formula = '`Litter spe. richness (#/dm2)`~sp*heterogeneity')
mod |> summary()

mod = lm(data = df, 
         formula = '`C decomposition rate (%)`~sp*heterogeneity')
mod |> summary()

mod = lm(data = df, 
         formula = '`C decomposition rate SD (%)`~sp*heterogeneity')
mod |> summar()

mod = lm(data = df, 
         formula = '`C loss (g/dm2)`~sp*heterogeneity')
mod |> summary()
