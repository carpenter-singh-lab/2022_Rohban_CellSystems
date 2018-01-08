rm(list = ls())

library(dplyr)
library(stringr)

load("../results/master/2017-12-12_b75637fc/workspace_MOApred.RData")

pthw <- gene.compound.cr.pr %>% 
  select(Var1, Name, MOA, Target) %>% 
  unique()
  
data.frame(compound = .compound, p.value = p.vals) %>%
  filter(!compound %in% compound.blacklist) %>%
  mutate(p.value = p.adjust(p.value, "BH")) %>%
  arrange(p.value) %>%
  mutate(p.value = round(p.value, 4)) %>%
  left_join(., pthw, by = c("compound" = "Var1")) %>%
  rename(adjusted.p.value = p.value) %>%
  select(Name, MOA, adjusted.p.value) %>%
  htmlTable::htmlTable()
