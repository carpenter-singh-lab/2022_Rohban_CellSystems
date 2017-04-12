rm(list = ls())
library(dplyr)
library(stringr)

load("../input/Gustafsdottir/Initial_analysis.RData")
moa <- read.csv("../input/moas.txt", sep = "\t")

cmpd.names <- data.frame(Name = 
                           Pf.full$Image_Metadata_SOURCE_COMPOUND_NAME %>%
                           unique) %>% unique
cmpd.names.lower <- cmpd.names
cmpd.names.lower$Name <- lapply(cmpd.names$Name, function(x) str_to_lower(x)) %>% unlist

cmpd.names.san <- cmpd.names.lower
cmpd.names.san$Name <- lapply(cmpd.names.san$Name, function(x) x %>% 
                           str_replace_all(., "-", "") %>% 
                           str_replace_all(., " ", "") %>% 
                           str_replace_all(., "\\+", "") %>% 
                           str_replace_all(., "\\(", "") %>% 
                           str_replace_all(., "\\)", "") %>% 
                           str_replace_all(., ",", "")) %>% unlist
cmpd.names.san <- cbind(cmpd.names.san, Name.orig = cmpd.names$Name)
  
moa.san <- moa
moa.san$Name <- lapply(moa.san$Name, function(x) x %>% 
                                str_replace_all(., "-", "") %>% 
                                str_replace_all(., " ", "") %>% 
                                str_replace_all(., "\\+", "") %>% 
                                str_replace_all(., "\\(", "") %>% 
                                str_replace_all(., "\\)", "") %>% 
                                str_replace_all(., ",", "")) %>% unlist

moa.san.join <- plyr::join(moa.san[,c("Name", "MOA", "Target")],
                           cmpd.names.san, by = "Name", type = "right") 
moa.san.join <- moa.san.join %>% dplyr::select(one_of(c("Name.orig", "MOA", "Target")))
colnames(moa.san.join) <- c("Name", "MOA", "Target")
moa.san.join <- moa.san.join %>% unique
moa.san.join %>% dplyr::filter(!is.na(MOA) & MOA != "") %>% NROW()
readr::write_csv(moa.san.join, "MOAs.csv")

