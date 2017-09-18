rm(list = ls())

library(dplyr)
library(stringr)

toxic.small <- "BRD-A49680073-001-01-9"
toxic.big <- "BRD-K51575138-001-01-7"
unq <- "BRD-A34806832-001-02-7"
thrd <- "BRD-K00259736-001-03-2"
frth <- "BRD-A18497530-001-04-6"
big.spec <- "BRD-K37798499-001-02-5"
btm <- "BRD-K41853443-001-01-1"
sbtm <- "BRD-K40742111-001-02-6"
adown <- "BRD-K96799727-001-01-7"
brd <- adown


Pf <- readRDS("../results/master/2017-09-05_da5f3073/Pf_bio_new_all.rds")

u <- Pf %>% 
  filter(Metadata_broad_sample == brd) %>%
  select(Metadata_Plate, Metadata_Well) %>%
  unique 

str <- sprintf("aws s3 sync s3://imaging-platform/projects/2015_Bray_GigaScience/CDRP/images/%s ~/work/projects/gene_compound/results/manual/toxic_phenotypes --exclude '*' --include '*_%s_s5_*'", u$Metadata_Plate[1], str_to_lower(u$Metadata_Well[1]))

print(str)
