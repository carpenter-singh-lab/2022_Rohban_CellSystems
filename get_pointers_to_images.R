rm(list = ls())
gc()

library(dplyr)
library(stringr)

brds <- c("BRD-K54330070-001-05-5")
pert.type <- "compound"
#brd <- adown
# find . -name "combined.png" -exec rm -rf {} \;

if (pert.type == "gene") {
  load("../input/TA/Initial_analysis_workspace_new.RData")  
  Pf <- Pf_org$data %>%
    mutate(Metadata_broad_sample = Treatment, 
           Metadata_Plate = Plate_, 
           Metadata_Well = Well) 
  
} else {
  Pf.1 <- readRDS("../results/master/2017-09-05_da5f3073/Pf_bio_new_all.rds")
  
  col.names <- colnames(Pf.1)
  col.names <- col.names[which(str_detect(col.names, "Metadata_"))]
  
  Pf.1 <- Pf.1 %>% select(one_of(col.names))
  
  Pf.2 <- readRDS("../results/master/2017-09-05_da5f3073/Pf_DOS_new_all.rds")
  Pf.2 <- Pf.2 %>% select(one_of(col.names))
  
  Pf <- rbind(Pf.1, Pf.2)
}

for (brd in brds) {

  u <- Pf %>% 
    filter(Metadata_broad_sample == brd) %>%
    select(Metadata_Plate, Metadata_Well) %>%
    unique
    
  if (pert.type == "gene")  {
    u <- u %>% filter(Metadata_Plate != "41749")
  }
  
  u <- u %>%
    slice(sample(1:NROW(u), 1))
  
  if (pert.type == "gene") {
    str <- sprintf("aws s3 sync s3://imaging-platform/projects/2011_07_13_TargetAccelerator_CancerProgram_MPG/SIGMA2_Pilot_2013_10_11/images/%s ~/work/projects/gene_compound/results/manual/misc/gene/%s --exclude '*' --include '*_%s_s5_*' --exclude '*_thumb*'", u$Metadata_Plate[1], brd, str_to_lower(u$Metadata_Well[1]))
  } else {
    str <- sprintf("aws s3 sync s3://imaging-platform/projects/2015_Bray_GigaScience/CDRP/images/%s ~/work/projects/gene_compound/results/manual/misc/compound/%s --exclude '*' --include '*_%s_s5_*'", u$Metadata_Plate[1], brd, str_to_lower(u$Metadata_Well[1]))
  }
  
  cat(str, sep = "\n")
  
}
