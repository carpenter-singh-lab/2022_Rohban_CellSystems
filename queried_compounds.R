library(dplyr)
library(stringr)
source("read_dataset.R")

Pf.repurp <- read.dataset("Repurposing", dose.closest = NULL)

a <- Pf.repurp$data %>%
  select(one_of(c(Pf.repurp$feat_cols, 
                  "Metadata_Treatment", 
                  "Metadata_pert_iname"))) %>% 
  group_by(Metadata_Treatment, 
           Metadata_pert_iname) %>% 
  summarise_each(funs("mean"))

cr <- cor(a[,Pf.repurp$feat_cols] %>% t)

colnames(cr) <- a$Metadata_Treatment
rownames(cr) <- a$Metadata_Treatment

cr.melt <- cr %>% 
  reshape2::melt() %>% 
  left_join(., 
            Pf.repurp$data[,c("Metadata_pert_iname", "Metadata_moa", "Metadata_Treatment")] %>% unique, 
            by = c("Var1" = "Metadata_Treatment")) %>%
  left_join(., 
            Pf.repurp$data[,c("Metadata_pert_iname", "Metadata_moa", "Metadata_Treatment")] %>% unique, 
            by = c("Var2" = "Metadata_Treatment")) %>%
  filter(Metadata_pert_iname.x != Metadata_pert_iname.y) %>%
  filter(Metadata_pert_iname.x %in% c("clofoctol",
                     "disulfiram",
                     "felbamate",
                     "KPT-330",
                     "MK-2206",
                     "mozavaptan",
                     "STA-5326",
                     "tosedostat")) 

#View(cr.melt)
u <- cr.melt %>% 
  filter(Metadata_pert_iname.y != "dmso") %>%
  arrange(-value) %>%
  group_by(Metadata_pert_iname.x) %>%
  slice(1:50) %>%
  ungroup %>%
  group_by(Metadata_pert_iname.x, Metadata_pert_iname.y, Metadata_moa.x, Metadata_moa.y) %>%
  tally() %>%
  ungroup() %>%
  arrange(-n) %>%
  group_by(Metadata_pert_iname.x) %>%
  slice(1:10) %>%
  ungroup()
  
u %>% htmlTable::htmlTable()
