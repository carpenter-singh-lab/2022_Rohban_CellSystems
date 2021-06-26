## rep corr. against cell density
rm(list = ls())

library(dplyr)
library(stringr)
library(ggplot2)
library(progress)

source("rep.corr.func.R")
  
Pf <- readRDS("../results/master/2017-04-20_425d653/Pf_bio_new.rds")
ft <- colnames(Pf)
meta.col <- ft[which(str_detect(ft, "Metadata_"))]
ft <- ft[which(!str_detect(ft, "Metadata_"))]

Pf.1 <- readRDS("../results/master/2017-09-05_da5f3073/Pf_bio_new_all.rds")
Pf.2 <- readRDS("../results/master/2017-09-05_da5f3073/Pf_DOS_new_all.rds")

Pf <- rbind(Pf.1, Pf.2)

Pf <- Pf %>% select(one_of(c(meta.col, ft)))

cell.count <- readr::read_csv("../results/master/2017-11-14_7a536be7/toxicity.csv")
colnames(cell.count) <- c("Metadata_Plate", "Metadata_Well", "total.dna")

dt <- Pf %>% 
  left_join(., cell.count, by = c("Metadata_Plate", "Metadata_Well")) %>%
  filter(Metadata_broad_sample == "DMSO") %>%
  select(total.dna) %>%
  as.matrix() %>%
  as.vector()

dt.grp <- Pf %>% 
  left_join(., cell.count, by = c("Metadata_Plate", "Metadata_Well")) %>%
  filter(Metadata_broad_sample == "DMSO") %>%
  group_by(Metadata_Plate) %>%
  summarise(total.dna.median = median(total.dna), total.dna.mad = mad(total.dna)) %>%
  ungroup()
  
mn <- mean(dt)
sdv <- sd(dt)

cell.counts <- Pf %>% 
  left_join(., cell.count, by = c("Metadata_Plate", "Metadata_Well")) 

pb <- progress_bar$new(total = NROW(dt.grp))

for (pl in dt.grp$Metadata_Plate) {
  mn <- dt.grp %>% filter(Metadata_Plate == pl) %>% select(total.dna.median) %>% as.matrix() %>% as.vector()
  sd <- dt.grp %>% filter(Metadata_Plate == pl) %>% select(total.dna.mad) %>% as.matrix() %>% as.vector()
  
  cell.counts <- cell.counts %>%
    mutate(total.dna = ifelse(Metadata_Plate == pl, (total.dna - mn)/sd, total.dna))
  
  pb$tick()
}

cell.counts <- cell.counts %>%
  group_by(Metadata_broad_sample) %>%
  summarise(total.dna = mean(total.dna, na.rm = T)) %>%
  ungroup() 

Pf <- Pf %>% 
  filter(Metadata_broad_sample != "DMSO") 

u <- rep.cor(list(data = Pf, feat_cols = ft), grp.var = "Metadata_broad_sample", feat.var = ft)
thr <- non.rep.cor(list(data = Pf, feat_cols = ft), grp.var = "Metadata_broad_sample", feat.var = ft, quant = 0.95)

print(length(which(u$cr > thr))/NROW(u))

u.joined <- u %>% 
  left_join(., cell.counts, by = "Metadata_broad_sample")

toxic.cmpd <- u.joined %>%
  filter(cr > thr & total.dna < -3)

Pf.collapsed <- Pf %>%
  filter(Metadata_broad_sample %in% toxic.cmpd$Metadata_broad_sample) %>%
  group_by(Metadata_broad_sample) %>%
  summarise_at(.funs = mean, .vars = ft) %>%
  ungroup() 

write.table(Pf.collapsed, "Pf_toxic_collapsed.txt", sep = "\t", row.names = F)
