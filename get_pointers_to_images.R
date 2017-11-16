rm(list = ls())
gc()

library(dplyr)
library(stringr)

brds <- c("BRD-K11018734-001-01-5", "BRD-K15475572-001-01-8", "BRD-A16504062-001-07-0",
          "BRD-K70510768-001-07-5", "BRD-K33368614-001-05-4", "BRD-K21360345-001-05-1",
          "BRD-K39396484-001-05-6", "BRD-K75464194-001-04-1", "BRD-K46123637-001-01-0",
          "BRD-K43438744-001-02-3", "BRD-A31650268-001-05-8", "BRD-K02764365-001-08-6",
          "BRD-K05789599-001-11-2", "BRD-K18235559-001-01-0", "BRD-K25646871-001-01-0", 
          "BRD-K06543683-001-03-6", "BRD-K64581877-001-01-5", "BRD-K40636143-001-08-9",
          "BRD-K08307026-001-01-4", "BRD-K68917053-001-08-8", "BRD-K71387353-001-01-1", 
          "BRD-K18350116-001-01-3", "BRD-K80045183-001-01-0", "BRD-K26053980-001-01-6", 
          "BRD-K90918886-001-01-6", "BRD-K79497053-001-01-3", "BRD-K29441085-001-01-4")

#brd <- adown

Pf.1 <- readRDS("../results/master/2017-09-05_da5f3073/Pf_bio_new_all.rds")

col.names <- colnames(Pf.1)
col.names <- col.names[which(str_detect(col.names, "Metadata_"))]

Pf.1 <- Pf.1 %>% select(one_of(col.names))

Pf.2 <- readRDS("../results/master/2017-09-05_da5f3073/Pf_DOS_new_all.rds")
Pf.2 <- Pf.2 %>% select(one_of(col.names))

Pf <- rbind(Pf.1, Pf.2)

for (brd in brds) {

  u <- Pf %>% 
    filter(Metadata_broad_sample == brd) %>%
    select(Metadata_Plate, Metadata_Well) %>%
    unique 
  
  str <- sprintf("aws s3 sync s3://imaging-platform/projects/2015_Bray_GigaScience/CDRP/images/%s ~/work/projects/gene_compound/results/manual/toxic_phenotypes/%s --exclude '*' --include '*_%s_s5_*'", u$Metadata_Plate[1], brd, str_to_lower(u$Metadata_Well[1]))
  
  cat(str, sep = "\n")
  
}
