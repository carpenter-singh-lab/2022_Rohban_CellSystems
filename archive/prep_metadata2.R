rm(list = ls())
library(dplyr)
library(stringr)

c <- readr::read_csv("../../ta/analysis/input/CDP2/cdp2.metadata.csv")
meta <- c %>% unique
meta1 <- meta
cl <- colnames(meta1)
cl.new <- lapply(cl, function(x) str_replace(x, "Image_Metadata_", "")) %>% unlist
colnames(meta1) <- cl.new

colnames(meta1)[4] <- "broad_sample"

meta1 <- meta1 %>% dplyr::mutate(mmoles_per_liter = CPD_MMOL_CONC)
meta1 <- meta1 %>% dplyr::mutate(solvent = "DMSO")
meta1 <- meta1 %>% dplyr::mutate(isControl = ifelse(ASSAY_WELL_ROLE == "mock", T, F))
meta1 <- meta1 %>% dplyr::mutate(ASSAY_WELL_ROLE = ifelse(isControl, "mock", "treated"))
meta1 <- meta1 %>% dplyr::mutate(well_position = str_to_lower(CPD_WELL_POSITION))
meta1 <- meta1 %>% dplyr::mutate(CPD_Plate_Map_Name = CPD_PLATE_MAP_NAME)
meta1 <- meta1 %>% dplyr::select(-isControl, -CPD_WELL_POSITION,
                                 -CPD_WELL_POSITION_1, -Wave, -Site, -CLSID, -UserStem,
                                 -Time, -CPD_PLATE_MAP_NAME, -CPD_MMOL_CONC, -PlateID_1,
                                 -PlateID, -CPD_WELL_POSITION)
meta1 <- meta1 %>% unique

v <- meta1$CPD_Plate_Map_Name %>% unique

for (pl in v) {
  meta1 %>% dplyr::filter(CPD_Plate_Map_Name == pl) %>% dplyr::select(-CPD_Plate_Map_Name) %>% unique %>%
    write.table(., sprintf("%s.txt", pl), sep = "\t", quote = T, row.names = F)
}
