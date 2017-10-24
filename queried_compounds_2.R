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

cmpds <- readr::read_csv("../input/flo_cmpds.csv", col_names = F)

selected.trt <- Pf.repurp$data %>%
  filter(str_to_lower(Metadata_pert_iname) %in% str_to_lower(cmpds$X2)) %>%
  select(Metadata_Treatment) %>% 
  as.matrix() %>% 
  as.vector()

cr <- cr[selected.trt, ]
cr.norm.gc <- apply(cr, 2, function(x) (ecdf(x)(x))) 
rownames(cr.norm.gc) <- rownames(cr)
colnames(cr.norm.gc) <- colnames(cr)

cr.melt <- cr %>% 
  reshape2::melt() %>% 
  filter(Var1 %in% selected.trt & Var2 %in% selected.trt) %>%
  left_join(., 
            Pf.repurp$data[,c("Metadata_pert_iname", "Metadata_moa", "Metadata_Treatment")] %>% unique, 
            by = c("Var1" = "Metadata_Treatment")) %>%
  left_join(., 
            Pf.repurp$data[,c("Metadata_pert_iname", "Metadata_moa", "Metadata_Treatment")] %>% unique, 
            by = c("Var2" = "Metadata_Treatment")) 

cr.melt %>% 
  filter(Metadata_pert_iname.y != "dmso") %>%
  group_by(Metadata_pert_iname.x, Metadata_pert_iname.y) %>%
  summarise(value = max(value)) %>%
  arrange(-value) %>% 
  filter(Metadata_pert_iname.x > Metadata_pert_iname.y) %>%
  filter(Metadata_pert_iname.x == "dopamine" | Metadata_pert_iname.y == "dopamine") %>%
  htmlTable::htmlTable()  
