rm(list = ls())
library(dplyr)

#source("read_dataset.R")
#x <- read.dataset("CDRP", just.bioactives = F)
x1 <- readRDS("../results/master/2017-09-05_da5f3073/Pf_bio_new_all.rds")
x2 <- readRDS("../results/master/2017-09-05_da5f3073/Pf_DOS_new_all.rds")

x <- rbind(x1, x2)

y <- x %>% select("Metadata_Plate", "Metadata_Well", "Metadata_ASSAY_WELL_ROLE", "Metadata_broad_sample")
x <- y

y <- readr::read_csv("../input/Cell_counts.csv")

y <- y %>% 
  left_join(., x, by = c("Image_Metadata_Plate" = "Metadata_Plate", "Image_Metadata_Well" = "Metadata_Well"))

readr::write_csv(y, "Cell_counts_annotated.csv")
