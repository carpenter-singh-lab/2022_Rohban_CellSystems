library(dplyr)
library(stringr)

cmpd.name <- "sb-202190"

moa2 <- readr::read_csv("../input/CDP2/MOA_annot2.csv")
moa2 <- moa2 %>%
  mutate(CPD_NAME = str_to_lower(CPD_NAME)) %>%
  mutate(CPD_NAME = ifelse(str_sub(CPD_NAME, 1, 4) == "brd-", Metadata_broad_sample, CPD_NAME))

moa2 <- moa2 %>%
  group_by(Name) %>%
  slice(1) %>%
  ungroup()

moa2 %>%
  filter(str_to_lower(Name) == cmpd.name) %>%
  select(Metadata_broad_sample) %>%
  unique() %>%
  as.matrix() %>%
  as.vector() %>%
  print