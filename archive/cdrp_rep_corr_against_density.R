## rep corr. against cell density

library(dplyr)
library(stringr)
library(ggplot2)

source("rep.corr.func.R")

Pf.1 <- readRDS("../results/master/2017-09-05_da5f3073/Pf_bio_new_all.rds")
Pf.2 <- readRDS("../results/master/2017-09-05_da5f3073/Pf_DOS_new_all.rds")

Pf <- rbind(Pf.1, Pf.2)

cell.count <- readr::read_csv("../results/master/2017-04-21_425d6537/Cell_counts.csv")
colnames(cell.count) <- c("Metadata_Plate", "Metadata_Well", "cell.count")

dt <- Pf %>%
  left_join(., cell.count, by = c("Metadata_Plate", "Metadata_Well")) %>%
  filter(Metadata_broad_sample == "DMSO") %>%
  select(cell.count) %>%
  as.matrix() %>%
  as.vector()

mn <- mean(dt, na.rm = T)
sdv <- sd(dt, na.rm = T)

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

u.summ <- u.joined %>%
  mutate(pos.cr = (cr > thr),
         pos.cell.count = (cell.count > 0))  %>%
  group_by(pos.cr, pos.cell.count) %>%
  tally() %>%
  ungroup() %>%
  mutate(freq = n / sum(n)) %>%
  select(-n)

g <- ggplot(u.joined, aes(x = cell.count, y = cr)) +
  geom_point(size = 0.5) + xlab("cell count z-score") +
  ylab("median replicate corr.") +
  geom_hline(yintercept = thr, color = "red") +
  geom_vline(xintercept = 0, color = "blue") +
  theme_bw() +
  geom_text(aes(2, 1, label = paste(as.character(round(u.summ[which(u.summ$pos.cr &
                                                                         u.summ$pos.cell.count) , "freq"] * 100,
                                                           0)), "%")), size = 6, color = "gray50") +
  geom_text(aes(-2, -0.40, label = paste(as.character(round(u.summ[which(!u.summ$pos.cr &
                                                                           !u.summ$pos.cell.count) , "freq"] * 100,
                                                           0)), "%")), size = 6, color = "gray50") +
  geom_text(aes(2, -0.40, label = paste(as.character(round(u.summ[which(!u.summ$pos.cr &
                                                                          u.summ$pos.cell.count) , "freq"] * 100,
                                                           0)), "%")), size = 6, color = "gray50") +
  geom_text(aes(-2, 1, label = paste(as.character(round(u.summ[which(u.summ$pos.cr &
                                                                          !u.summ$pos.cell.count) , "freq"] * 100,
                                                             0)), "%")), size = 6, color = "gray50")

ggsave("rep_cor_vs_density.png", g)
