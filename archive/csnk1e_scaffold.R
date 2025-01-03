rm(list = ls())
library(dplyr)
library(stringr)

cr.melt <- readRDS("../results/master/2017-10-11_7d18c89d/cr_melt_cp_canon.rds")
cr.melt.l1000 <- readRDS("../results/master/2017-10-11_65efa2d7/cr_melt_l1000.rds")
isomers <- read.csv("../input/csnk1e_scaffold")

cr.melt %>%
  mutate(Var1.trm = Var1) %>%
  filter(Var1.trm %in% isomers$BroadID) %>%
  arrange(value) %>%
  group_by(Var1.trm) %>%
  summarise(match.to.CSNK1E = (Var2[1] == "CSNK1E_WT.2"),
            score = value[1],
            specificity = (value[1] - value[2])/value[2],
            rank = which(Var2 == "CSNK1E_WT.2")) %>%
  ungroup() %>%
  mutate(specificity = ifelse(rank == 1 & -score >= 0.4, specificity, 0)) %>%
  arrange(-specificity) %>%
  left_join(., cr.melt.l1000 %>% filter(Var2 == "CSNK1E_WT.2"),
            by = c("Var1.trm" = "Var1")) %>%
  rename(BroadID = Var1.trm,
         CellPainting_score = score,
         L1000_score = value,
         CellPainting_specificity = specificity) %>%
  select(BroadID, CellPainting_score, CellPainting_specificity, L1000_score) %>%
  readr::write_csv(., "csnk1e_scaffolds_matches.csv")


