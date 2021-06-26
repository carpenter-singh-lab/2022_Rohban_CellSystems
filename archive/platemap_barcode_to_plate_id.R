rm(list = ls())
library(dplyr)
library(stringr)

c <- readRDS("../input/CDP2/Pf.rds")
d <- c$data[,c("Image_Metadata_CPD_Plate_Map_Name", "Image_Metadata_Plate")]

d <- unique(d)
colnames(d) <- c("Plate_Map_Name", "Assay_Plate_Barcode")

readr::write_csv(d, "barcode_platemap.csv")
