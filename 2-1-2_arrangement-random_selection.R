# R version 4.3.1 (2023-06-16)
# Platform: aarch64-apple-darwin20 (64-bit)
# Running under: macOS Sonoma 14.6.1

# Packages
# tidyverse_2.0.0
# arrangements_1.1.9
rm(list = ls())
set.seed(1603)

library(tidyverse)
library(arrangements)

# Counting number of potential combinations and arrangements
# 1 species out of 8
ncombinations(8, 1)
arrangement.sp.1 = combinations(8,1)

# 2 species out of 8
ncombinations(8, 2)
arrangement.sp.2 = combinations(8, 2)

# 4 species out of 8
ncombinations(8, 4)
npermutations(8, 4)
a = permutations(8, 4) %>%
  data.frame() %>%
  mutate(composition = 2^X1 + 2^X2 + 2^X3 + 2^X4) %>%
  mutate(arrangement = 1:nrow(.))
  
a.select = c()
for(i in unique(a$composition)){
  a.select = 
    c(
      a.select,
      a %>% 
        filter(composition == i) %>%
        sample_n(15) %>% 
        select(arrangement)
    )
}

# Selecting a subset of 1000 permutations
a.select = 
  a.select %>% 
  unlist %>%
  sample(size = 1000) %>%
  as.numeric()

arrangement.sp.4 = 
  a %>% 
  filter(arrangement %in% a.select) 

arrangement.sp.4 %>%
  group_by(composition) %>%
  summarise(n = n()) %>% 
  arrange(n) %>%
  ggplot(aes(x = n)) + 
  geom_density() + 
  scale_x_continuous(breaks = 0:20)

# 8 species out of 8
npermutations(8,8)
arrangement.sp.8 = 
  permutations(8, 8) %>%
  data.frame() %>%
  mutate(composition = 2^X1 + 2^X2 + 2^X3 + 2^X4) %>%
  mutate(arrangement = 1:nrow(.)) %>%
  # Selecting a subset of 1000 permutations
  sample_n(1000)

# Saving the selected arrangements
arrangements = list(arrangement.sp.1, arrangement.sp.2, arrangement.sp.4, arrangement.sp.8)
write_rds(arrangements, file = '2-1-2_forests-arrangements.RDS')
#### END ####