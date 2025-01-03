rm(list = ls())
gc()

library(dplyr)
library(stringr)
source("read_dataset.R")

brds <- c("BRD-K69776681-001-03-8")
pert.type <- "compound"
ds <- "Repurp"
#brd <- adown
# find . -name "combined.png" -exec rm -rf {} \;

if (pert.type == "gene") {
  load("../input/TA/Initial_analysis_workspace_new.RData")
  Pf <- Pf_org$data %>%
    mutate(Metadata_broad_sample = Treatment,
           Metadata_Plate = Plate_,
           Metadata_Well = Well)

} else if (ds == "CDRP") {
  Pf.1 <- readRDS("../results/master/2017-09-05_da5f3073/Pf_bio_new_all.rds")

  col.names <- colnames(Pf.1)
  col.names <- col.names[which(str_detect(col.names, "Metadata_"))]

  Pf.1 <- Pf.1 %>% select(one_of(col.names))

  Pf.2 <- readRDS("../results/master/2017-09-05_da5f3073/Pf_DOS_new_all.rds")
  Pf.2 <- Pf.2 %>% select(one_of(col.names))

  Pf <- rbind(Pf.1, Pf.2)
} else {
  Pf <- read.dataset("Repurposing")
  Pf <- Pf$data
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
  } else if (ds == "CDRP") {
    str <- sprintf("aws s3 sync s3://imaging-platform/projects/2015_Bray_GigaScience/CDRP/images/%s ~/work/projects/gene_compound/results/manual/misc/compound/%s --exclude '*' --include '*_%s_s5_*'", u$Metadata_Plate[1], brd, str_to_lower(u$Metadata_Well[1]))
  } else {
    str1 <- system("aws s3 ls s3://imaging-platform/projects/2015_10_05_DrugRepurposing_AravindSubramanian_GolubLab_Broad/2016_04_01_a549_48hr_batch1/images/", intern = T)

    str2 <- str1[which(stringr::str_detect(str1, u$Metadata_Plate[1]))]
    str3 <- stringr::str_split(str2, "PRE ")[[1]][2]
    str3 <- paste0(str3, "Images")

    i1 <- which(letters == stringr::str_sub(u$Metadata_Well[1], 1, 1))
    i1.s <- i1 %>% as.character()
    if (i1 < 10) {
      i1.s <- paste0("0", i1.s)
    }
    i2 <- stringr::str_sub(u$Metadata_Well[1], 2, 3)

    str <- sprintf("aws s3 sync s3://imaging-platform/projects/2015_10_05_DrugRepurposing_AravindSubramanian_GolubLab_Broad/2016_04_01_a549_48hr_batch1/images/%s ~/work/projects/gene_compound/results/manual/misc/compound/%s --exclude '*' --include 'r%sc%sf05*'", str3, brd, i1.s, i2)
  }

  cat(str, sep = "\n")

}
