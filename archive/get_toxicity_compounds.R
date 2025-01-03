## prepare the toxicity index

library(dplyr)
library(stringr)
library(foreach)

lst <- list.dirs("../input/toxicity")
toxicity <- foreach (dr = lst[2:length(lst)]) %do% {
  lst.sub <- list.files(dr)

  readr::read_csv(paste0(dr, "/", lst.sub[1])) %>%
    select(-sum.y, -sum) %>%
    rename(Metadata_Plate = Image_Metadata_Plate, Metadata_Well = Image_Metadata_Well)
}

toxicity <- do.call(rbind, toxicity)

readr::write_csv(toxicity, "toxicity.csv")
