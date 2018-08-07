rm(list = ls())
gc()

library(dplyr)
library(stringr)

source("read_dataset.R")

a <- read.dataset("CDRP", just.bioactives = F)
a$data <- a$data %>% 
  mutate(Metadata_broad_sample = str_sub(Metadata_broad_sample, 1, 13))
  
cell.count <- read.csv("../input/Cell_counts.csv") %>%
  rename(Metadata_Plate = Image_Metadata_Plate, Metadata_Well = Image_Metadata_Well) %>%
  left_join(a$data %>% 
              select(Metadata_Plate, Metadata_Well, Metadata_broad_sample)) %>%
  mutate(Metadata_broad_sample = str_sub(Metadata_broad_sample, 1, 13)) 

v <- cell.count %>%
  filter(Metadata_broad_sample == "DMSO") %>%
  select(Image_Count_Cells) %>% 
  as.matrix() %>% 
  as.vector()

v.mean <- mean(v)
v.std <- sd(v)

cell.count <- cell.count %>%
  group_by(Metadata_broad_sample) %>%
  summarise(cell.count = (mean(Image_Count_Cells) - v.mean)/v.std)

b <- readRDS("../results/master/2017-10-11_7d18c89d/cr_melt_cp_canon.rds") %>%
  mutate(Var1 = str_sub(Var1, 1, 13))
  
brd.list <- read.csv("../input/ppargc1_jon_followup.csv", header = F) %>% mutate(V1 = str_sub(V1, 1, 13)) %>% as.matrix() %>% as.vector()

a$data %>%
  filter(Metadata_broad_sample %in% brd.list) %>%
  select(Metadata_broad_sample, Cells_AreaShape_Area, Cells_Intensity_MeanIntensity_AGP, Cells_Intensity_StdIntensityEdge_Mito) %>%
  group_by(Metadata_broad_sample) %>%
  summarise(Cells_AreaShape_Area = median(Cells_AreaShape_Area), 
            Cells_Intensity_MeanIntensity_AGP = median(Cells_Intensity_MeanIntensity_AGP), 
            Cells_Intensity_StdIntensityEdge_Mito = median(Cells_Intensity_StdIntensityEdge_Mito)) %>%
  left_join(cell.count, by = "Metadata_broad_sample") %>% 
  left_join(b %>% filter(Var2 == "PPARGC1A_WT.2"), by = c("Metadata_broad_sample" = "Var1")) %>%
  left_join(a$data %>% select(Metadata_broad_sample, Metadata_pert_iname) %>% unique, by = "Metadata_broad_sample") %>%
  mutate(Metadata_pert_iname = ifelse(str_sub(Metadata_pert_iname, 1, 4) == "brd-", "", Metadata_pert_iname)) %>%
  select(-Var2) %>%
  mutate(Cells_AreaShape_Area = round(Cells_AreaShape_Area, 2),
         Cells_Intensity_MeanIntensity_AGP = round(Cells_Intensity_MeanIntensity_AGP, 2), 
         Cells_Intensity_StdIntensityEdge_Mito = round(Cells_Intensity_StdIntensityEdge_Mito, 2),
         cell.count = round(cell.count, 2), 
         value = round(value, 2)) %>%
  arrange(value) %>%
  htmlTable::htmlTable()

