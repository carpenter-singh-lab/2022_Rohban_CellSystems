## rep corr. against cell density
rm(list = ls())

library(dplyr)
library(stringr)
library(ggplot2)

source("rep.corr.func.R")

Pf <- readRDS("../results/master/2017-09-05_da5f3073/Pf_bio_new_all.rds")

cell.count <- readr::read_csv("../results/master/2017-04-21_425d6537/Cell_counts.csv")
colnames(cell.count) <- c("Metadata_Plate", "Metadata_Well", "cell.count")

dt <- Pf %>% 
  left_join(., cell.count, by = c("Metadata_Plate", "Metadata_Well")) %>%
  filter(Metadata_broad_sample == "DMSO") %>%
  select(cell.count) %>%
  as.matrix() %>%
  as.vector()

mn <- mean(dt)
sdv <- sd(dt)

cell.counts <- Pf %>% 
  left_join(., cell.count, by = c("Metadata_Plate", "Metadata_Well")) %>%
  group_by(Metadata_broad_sample) %>%
  summarise(cell.count = (mean(cell.count, na.rm = T) - mn)/sdv) %>%
  ungroup() 

Pf <- Pf %>% 
  filter(Metadata_broad_sample != "DMSO")

ft <- colnames(Pf)
ft <- ft[which(!str_detect(ft, "Metadata_"))]
u <- rep.cor(list(data = Pf, feat_cols = ft), grp.var = "Metadata_broad_sample", feat.var = ft)
thr <- non.rep.cor(list(data = Pf, feat_cols = ft), grp.var = "Metadata_broad_sample", feat.var = ft, quant = 0.95)

print(length(which(u$cr > thr))/NROW(u))

u.joined <- u %>% 
  left_join(., cell.counts, by = "Metadata_broad_sample")

toxic.cmpd <- u.joined %>%
  filter(cr > thr & cell.count < -1)

Pf.collapsed <- Pf %>%
  filter(Metadata_broad_sample %in% toxic.cmpd$Metadata_broad_sample) %>%
  group_by(Metadata_broad_sample) %>%
  summarise_at(.funs = mean, .vars = ft)

readr::write_csv(Pf.collapsed, "Pf_toxic_collapsed.csv")
