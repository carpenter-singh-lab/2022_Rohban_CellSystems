rm(list = ls())
library(dplyr)
library(readr)
library(stringr)
library(foreach)
library(htmlTable)

corr <- readRDS("../results/master/2017-10-11_7d18c89d/cr_melt_cp_canon.rds")
metadata <- readr::read_csv("../input/CDP2/metadata_CDRP.csv")

corr <- corr %>% left_join(metadata, by = c("Var1" = "Metadata_broad_sample"))

corr.sub <- corr %>%
  filter(str_detect(Var2, "YAP1") | str_detect(Var2, "WWTR1")) %>%
  mutate(Metadata_pert_iname = ifelse(is.na(Metadata_pert_iname), Metadata_pert_iname2, Metadata_pert_iname)) %>%
  group_by(Var1, Metadata_moa, Metadata_pert_iname, Metadata_target) %>%
  summarise(value = median(value)) %>%
  ungroup()
  

pos.bio <- corr.sub %>% 
  filter(!is.na(Metadata_pert_iname) | !is.na(Metadata_moa) | !is.na(Metadata_target)) %>%
  arrange(-value) %>%
  slice(1:15) %>%
  filter(abs(value) > 0.4)

neg.bio <- corr.sub %>% 
  filter(!is.na(Metadata_pert_iname) | !is.na(Metadata_moa) | !is.na(Metadata_target)) %>%
  arrange(value) %>%
  slice(1:15) %>%
  filter(abs(value) > 0.4)

pos.dos <- corr.sub %>% 
  filter(!(!is.na(Metadata_pert_iname) | !is.na(Metadata_moa) | !is.na(Metadata_target))) %>%
  arrange(-value) %>%
  slice(1:30) %>%
  filter(abs(value) > 0.4)

neg.dos <- corr.sub %>% 
  filter(!(!is.na(Metadata_pert_iname) | !is.na(Metadata_moa) | !is.na(Metadata_target))) %>%
  arrange(value) %>%
  slice(1:30) %>%
  filter(abs(value) > 0.4)

all.cmpds <- c(pos.bio %>% select(Var1) %>% as.matrix() %>% as.vector(), 
  neg.bio %>% select(Var1) %>% as.matrix() %>% as.vector(),
  pos.dos %>% select(Var1) %>% as.matrix() %>% as.vector(),
  neg.dos %>% select(Var1) %>% as.matrix() %>% as.vector()) %>%
  unique()

all.cmpds %>% length()
all.cmpds %>% cat(., sep = "\n")

avai <- readr::read_csv("../results/manual/CDRP/yap/available_samples.csv", col_names = F)
avai <- avai %>% as.matrix() %>% as.vector()
avai <- avai[seq(from = 1, to = length(avai), by = 3)]

cell.count <- readr::read_csv("../input/CDP2/Cell_counts.csv")

plate.barcode <- readr::read_csv("../input/CDP2/barcode_platemap.csv")

read.all <- function(path) {
  lst <- list.files(path)  
  meta <- foreach (lsti = lst, .combine = rbind) %do% { 
    path.ext <- paste0(path, "/", lsti)
    pl.map <- str_split(lsti, "\\.txt")[[1]][1]
    read.csv(path.ext, sep = "\t") %>% 
      mutate(plate.map.name = pl.map)
  }
  return(meta)
}

plate.maps <- read.all("../input/CDP2/platemap")
plate.maps <- plate.maps %>%
  left_join(plate.barcode, by = c("plate.map.name" = "Plate_Map_Name")) 

plate.maps <- plate.maps %>%
  select(broad_sample, well_position, Assay_Plate_Barcode, ASSAY_WELL_ROLE) %>%
  rename(Metadata_broad_sample = broad_sample, Metadata_Well = well_position, Metadata_Plate = Assay_Plate_Barcode)

cell.count <- cell.count %>%
  left_join(plate.maps, by = c("Image_Metadata_Well" = "Metadata_Well", "Image_Metadata_Plate" = "Metadata_Plate"))

dmso.count <- cell.count %>% 
  filter(ASSAY_WELL_ROLE == "mock") %>%
  select(Image_Count_Cells) %>%
  as.matrix() %>%
  as.vector()

mn.dmso <- mean(dmso.count, na.rm = T)  
sd.dmso <- sd(dmso.count, na.rm = T)  

cell.count <- cell.count %>%
  mutate(Image_Count_Cells = (Image_Count_Cells - mn.dmso)/sd.dmso)

cell.count <- cell.count %>%
  filter(ASSAY_WELL_ROLE != "mock") %>%
  select(Metadata_broad_sample, Image_Count_Cells) %>%
  group_by(Metadata_broad_sample) %>%
  summarise(Image_Count_Cells = mean(Image_Count_Cells)) %>%
  ungroup()
  

corr.sub2 <- corr %>%
  filter(Var2 == "TRAF2_WT") %>%
  filter(Var1 %in% avai) %>%
  select(Var1, value)

corr.sub %>%
  filter(Var1 %in% avai) %>% 
  left_join(cell.count, by = c("Var1" = "Metadata_broad_sample")) %>%
  arrange(value) %>%
  left_join(corr.sub2, by = c("Var1")) %>%
  rename(corr.to.yap = value.x, corr.to.traf2 = value.y) %>%
  htmlTable::htmlTable()
