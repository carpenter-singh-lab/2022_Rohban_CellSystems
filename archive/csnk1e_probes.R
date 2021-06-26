rm(list = ls())
library(dplyr)
library(stringr)

cr.melt <- readRDS("../results/master/2017-10-11_7d18c89d/cr_melt_cp_canon.rds")
isomers <- read.csv("../input/csnk1e_stereoisomer.csv", header = F)

cr.melt %>% 
  mutate(Var1.trm = str_sub(Var1, 1, 13)) %>%
  filter(Var1.trm %in% isomers$V1) %>%
  arrange(value) %>%
  group_by(Var1.trm) %>%
  slice(1:2) %>%
  summarise(match.to.CSNK1E = (Var2[1] == "CSNK1E_WT.2"), 
            score = value[1], 
            specificity = (value[1] - value[2])/value[2]) %>%
  ungroup() 
