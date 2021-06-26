library(dplyr)
library(stringr)

source("rep.corr.func.R")

Pf <- readRDS("../results/master/2017-09-05_e2847e98/Pf_DOS_new_all.rds")

Pf <- Pf %>% 
  filter(Metadata_broad_sample != "DMSO")

ft <- colnames(Pf)
ft <- ft[which(!str_detect(ft, "Metadata_"))]
u <- rep.cor(list(data = Pf, feat_cols = ft), grp.var = "Metadata_broad_sample", feat.var = ft)
thr <- non.rep.cor(list(data = Pf, feat_cols = ft), grp.var = "Metadata_broad_sample", feat.var = ft, quant = 0.95)

print(length(which(u$cr > thr))/NROW(u))
