# R version 4.3.1 (2023-06-16)
# Platform: aarch64-apple-darwin20 (64-bit)
# Running under: macOS Sonoma 14.6.1
rm(list = ls())
set.seed(1603)

#### > 1. Packages and functions ####
# tidyverse_2.0.0

library(tidyverse)
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

#### > 1. Loading data ####
# Data handling
forests = readRDS(file = '2-1-3_species-randomization-8sp.RDS')
list.records = 
  list.files('simul/8/')
list.records = list.records[grep('output', list.records)]
x = list.records[1]
df = map_df(.x = list.records[-1], 
            function(x){
              data.frame(
                rep = x, 
                read.csv(paste0('simul/8/',x)))
            }) 
name.short = colnames(df)
colnames(df) = c(
  'Replicate','Design' , 'Total tree biomass (T)', 'Tree biomass SD (T)', 
  'Total litterfall (g)', 'Litterfall (g/dm2)', 'Litterfall SD (g/dm2)',
  'Litter spe. richness (#/dm2)','Litter species richness SD (#/dm2)',
  'C decomposition rate (%)', 'C decomposition rate SD (%)',
  'Mean nitrogen decomposition rate (%)', 'Nitrogen decomposition rate SD (%)',
  'C loss (g)','Total N loss (g)'
)

df$Replicate = 
  df$Replicate |> 
  str_remove_all('output-8sp_rep-') |> 
  str_remove_all('.csv') |> 
  as.numeric()

df.1 = df |>
  pivot_longer(cols = 3:15)

df.1$name = df.1$name |>
  factor(levels = c(
    'Total tree biomass (T)', 'Tree biomass SD (T)',
    'Total litterfall (g)', 'Litterfall (g/dm2)', 'Litterfall SD (g/dm2)',
    'Litter spe. richness (#/dm2)','Litter species richness SD (#/dm2)',
    'C decomposition rate (%)', 'C decomposition rate SD (%)',
    'Mean nitrogen decomposition rate (%)', 'Nitrogen decomposition rate SD (%)',
    'C loss (g)','Total N loss (g)'
  ))

# Calculate heterogeneity level
df.1$heterogeneity = NA
for(i in 1:length(forests$design)){
  df.1$heterogeneity[df.1$Design == i] = -(to.optimise(forests$design[[i]])+1)
}

df.1$value[df.1$name == 'C loss (g)'] = 
  df.1$value[df.1$name == 'C loss (g)']/(17*17*100)
df.2$name[df.2$name == 'C loss (g)'] = 'C loss (g/dm2)'
df.2 = df.1 |> filter(name %in% c('Total tree biomass (T)', 
                                  'Tree biomass SD (T)',
                                  'Litterfall (g/dm2)', 
                                  'Litterfall SD (g/dm2)',
                                  'Litter spe. richness (#/dm2)',
                                  'C decomposition rate (%)', 
                                  'C decomposition rate SD (%)',
                                  'C loss (g/dm2)'))
df.2$name = as.character(df.2$name)
df.2$name = factor(df.2$name, levels = c('Total tree biomass (T)', 
                                         'Tree biomass SD (T)',
                                         'Litterfall (g/dm2)', 
                                         'Litterfall SD (g/dm2)',
                                         'Litter spe. richness (#/dm2)',
                                         'C decomposition rate (%)', 
                                         'C decomposition rate SD (%)',
                                         'C loss (g/dm2)'))

df.3 = df.2 |> 
  filter(name %in% c('Total tree biomass (T)', 
                     'Tree biomass SD (T)'))

#### > 2. Plotting results ####
p.biomass = 
  ggplot(data = df.3, 
         aes(x = heterogeneity, 
             y = value, 
             color = factor(Design))) + 
  geom_jitter(data = df.3, 
              aes(x = heterogeneity, 
                  y = value, 
                  color = factor(Design)),
              size = .5, 
              alpha = .5, 
              width = 10) +
  geom_smooth(data = df.3|> 
                filter(!(Design %in% 2:5)), 
              aes(x = heterogeneity, 
                  y = value, 
                  color = factor(Design)),
              color = 'black', 
              method="loess", 
              level=0.90,
              size = .5) + 
  facet_wrap(vars(name), scales = 'free',ncol =  3) + 
  scale_color_manual(
    values = c(
      '1' = 'darkred',
      '2' = 'red', 
      '3' = 'darkblue',
      '4' = 'blue', '5' = 'gray70',
      '6' = 'gray60','7' = 'gray50','8' = 'gray40',
      '9' = 'gray30','10' = 'gray20','11' = 'gray10', '12' = 'black')) + 
  # scale_y_continuous(labels = scales::label_number_si()) +
  scale_x_continuous(breaks = c(-745, -670, -595,-380, 0),
                     labels = c("Blocks", 'Mini blocks',
                                "2 lines", 'Lines', 'Random')) +
  labs(x = "", 
       y = "Value", 
       title = "A. Tree biomass") +
  theme_bw() + 
  theme(legend.position = 'none', 
        axis.text.x = element_text(angle = 45, hjust = 1))
p.biomass 


df.3 = df.2 |> 
  filter(name %in% c('Litterfall (g/dm2)', 
                     'Litterfall SD (g/dm2)',
                     'Litter spe. richness (#/dm2)'))
p.litter = 
  ggplot(data = df.3, 
         aes(x = heterogeneity, 
             y = value, 
             color = factor(Design))) + 
  geom_jitter(data = df.3, 
                     aes(x = heterogeneity, 
                         y = value, 
                         color = factor(Design)),
                     size = .5, 
                     alpha = .5, 
              width = 10) +
  geom_smooth(data = df.3|> 
                filter(!(Design %in% 2:5)), 
              aes(x = heterogeneity, 
                  y = value, 
                  color = factor(Design)),
              color = 'black', 
              method="loess", 
              level=0.90,
              size = .5) + 
  facet_wrap(vars(name), scales = 'free',ncol =  3) + 
  scale_color_manual(
    values = c(
      '1' = 'darkred',
      '2' = 'red', 
      '3' = 'darkblue',
      '4' = 'blue', '5' = 'gray70',
      '6' = 'gray60','7' = 'gray50','8' = 'gray40',
      '9' = 'gray30','10' = 'gray20','11' = 'gray10', '12' = 'black')) + 
  # scale_y_continuous(labels = scales::label_number_si()) +
  scale_x_continuous(breaks = c(-745, -670, -595,-380, 0),
                     labels = c("Blocks", 'Mini blocks',
                                "2 lines", 'Lines', 'Random')) +
  labs(x = "", 
       y = "Value", 
       title = "B. Litterfall") +
  theme_bw() + 
  theme(legend.position = 'none', 
        axis.text.x = element_text(angle = 45, hjust = 1))
p.litter 

df.4 = df.2 |> filter(name %in% c('C decomposition rate (%)', 
                                  'C decomposition rate SD (%)',
                                  'C loss (g/dm2)'))

p.decomp = 
  ggplot(data = df.4, 
         aes(x = heterogeneity, 
             y = value, 
             color = factor(Design))) + 
  geom_jitter(data = df.4, 
              aes(x = heterogeneity, 
                  y = value, 
                  color = factor(Design)),
              size = .5, 
              alpha = .5,
              width = 10) +
  geom_smooth(data = df.4 |> 
                filter(!(Design %in% 2:4)), 
              aes(x = heterogeneity, 
                  y = value, 
                  color = factor(Design)),
              color = 'black', 
              method="loess", 
              level=0.90,
              size = .5) + 
  facet_wrap(vars(name), scales = 'free',ncol =  3) + 
  scale_color_manual(
    values = c(
      '1' = 'darkred',
      '2' = 'red', 
      '3' = 'darkblue',
      '4' = 'blue', '5' = 'gray70',
      '6' = 'gray60','7' = 'gray50','8' = 'gray40',
      '9' = 'gray30','10' = 'gray20','11' = 'gray10', '12' = 'black')) + 
  # scale_y_continuous(labels = scales::label_number_si()) +
  scale_x_continuous(breaks = c(-745, -670, -595,-380, 0),
                     labels = c("Blocks", 'Mini blocks',
                                "2 lines", 'Lines', 'Random')) +
  labs(x = "Tree species spatial heterogeneity", 
       y = "Value", 
       title = "C. Decomposition") +
  theme_bw() + 
  theme(legend.position = 'none', 
        axis.text.x = element_text(angle = 45, hjust = 1))
p.decomp  


Fig.2 = 
  ggpubr::ggarrange(
    ggpubr::ggarrange(p.biomass, 
                      ggplot(data=NULL) + theme_void(), 
                      nrow = 1, widths = c(2,1)),
    p.litter,
    p.decomp, 
    nrow = 3,
    align = 'hv')
Fig.2

#### 3. Saving results ####
ggsave(plot = Fig.2, 
       filename = '2-8-simulations-V2-8sp-based/3-1_Fig2.png',
       height = 20, width = 16, units = 'cm')

ggsave(plot = Fig.2, 
       filename = '2-8-simulations-V2-8sp-based/3-1_Fig2.svg',
       height = 25, width = 20, units = 'cm')
#### END ####