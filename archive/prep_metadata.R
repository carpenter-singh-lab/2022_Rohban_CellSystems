rm(list = ls())
library(dplyr)
library(stringr)

c <- readRDS("../input/CDP2/Pf.rds")
meta <- c$data[,c$factor_cols] %>% unique
meta1 <- meta %>% dplyr::filter(str_detect(Image_Metadata_CPD_Plate_Map_Name, "BIO"))
cl <- colnames(meta1)
cl.new <- lapply(cl, function(x) str_replace(x, "Image_Metadata_", "")) %>% unlist
colnames(meta1) <- cl.new

colnames(meta1)[1] <- "broad_sample"

meta1 <- meta1 %>% dplyr::mutate(mmoles_per_liter = 5)
meta1 <- meta1 %>% dplyr::mutate(solvent = "DMSO")
meta1 <- meta1 %>% dplyr::mutate(ASSAY_WELL_ROLE = ifelse(isControl, "mock", "treated"))
meta1 <- meta1 %>% dplyr::mutate(well_position = str_to_lower(Well))
meta1 <- meta1 %>% dplyr::select(-isControl, -WellUniqueID)
meta1 <- meta1 %>% dplyr::select(-Plate, -Well) %>% unique

v <- meta1$CPD_Plate_Map_Name %>% unique

for (pl in v) {
    meta1 %>% dplyr::filter(CPD_Plate_Map_Name == pl) %>% dplyr::select(-CPD_Plate_Map_Name) %>% unique %>%
    write.table(., sprintf("%s.txt", pl), sep = "\t", quote = T, row.names = F)
}
