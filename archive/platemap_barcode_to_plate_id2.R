rm(list = ls())
library(dplyr)
library(stringr)

c <- readr::read_csv("../../ta/analysis/input/CDP2/cdp2.metadata.csv")
meta <- c %>% unique

d <- c[,c("Image_Metadata_CPD_PLATE_MAP_NAME", "Image_Metadata_PlateID")]

d <- unique(d)
colnames(d) <- c("Plate_Map_Name", "Assay_Plate_Barcode")

readr::write_csv(d, "barcode_platemap.csv")
