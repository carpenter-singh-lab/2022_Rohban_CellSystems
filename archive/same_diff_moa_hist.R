rm(list = ls())
library(dplyr)
library(stringr)
library(doParallel)
source("moa_evaluations.R")

doParallel::registerDoParallel(cores = 4)

whiten <- F
top.p <- 0.01

cdrp <- readRDS("../results/master/2017-04-20_425d653/Pf_bio_new.rds")
feats <- colnames(cdrp)
feats <- feats[which(!str_detect(feats, "Metadata_"))]

if (whiten) {
  cdrp.dmso <- readRDS("../results/master/2017-09-05_da5f3073/Pf_bio_new_all.rds")
  cdrp.dmso <- cdrp.dmso %>% filter(Metadata_ASSAY_WELL_ROLE == "mock") %>% select(one_of(feats))
  mn <- apply(cdrp.dmso, 2, mean)
  cv <- cov(cdrp.dmso)
  ev <- eigen(cv)
  eps <- 10^-1
  W <- diag((ev$values + eps)^-0.5) %*% t(ev$vectors)
  cdrp[, feats] <- t(W %*% (apply(cdrp[, feats], 1, function(x) (x - mn))))
}

brd.to.name <- readr::read_csv("../input/CDP2/cdrp.cpd.meta.csv")
moa <- read.csv("../input/moas.txt", sep = "\t")

moa <- moa %>% 
  mutate(Name.cano = str_to_lower(Name)) %>%
  mutate(Name.cano = str_replace_all(Name.cano, "-", "")) %>%
  mutate(Name.cano = str_replace_all(Name.cano, " ", "")) %>%
  mutate(Name.cano = str_replace_all(Name.cano, "\\[", "")) %>%
  mutate(Name.cano = str_replace_all(Name.cano, "\\]", "")) %>%
  mutate(Name.cano = str_replace_all(Name.cano, "\\(", "")) %>%
  mutate(Name.cano = str_replace_all(Name.cano, "\\)", "")) %>%
  mutate(Name.cano = str_replace_all(Name.cano, "\\+", "")) 

cdrp.profiles <- cdrp %>% 
  group_by(Metadata_broad_sample) %>%
  summarise_at(.vars = feats, .funs = mean) %>%
  mutate(Metadata_broad_id_trunc = str_sub(Metadata_broad_sample, 1, 13)) %>%
  left_join(., brd.to.name, by = c("Metadata_broad_id_trunc" = "BROAD_CPD_ID")) %>% 
  filter(CPD_NAME != "")

cdrp.profiles <- cdrp.profiles %>% 
  mutate(Name.cano = str_to_lower(CPD_NAME)) %>%
  mutate(Name.cano = str_replace_all(Name.cano, "-", "")) %>%
  mutate(Name.cano = str_replace_all(Name.cano, " ", "")) %>%
  mutate(Name.cano = str_replace_all(Name.cano, "\\[", "")) %>%
  mutate(Name.cano = str_replace_all(Name.cano, "\\]", "")) %>%
  mutate(Name.cano = str_replace_all(Name.cano, "\\(", "")) %>%
  mutate(Name.cano = str_replace_all(Name.cano, "\\)", "")) %>%
  mutate(Name.cano = str_replace_all(Name.cano, "\\+", "")) 

moa2 <- readr::read_csv("../input/CDP2/MOA_annot2.csv")
moa2 <- moa2 %>% mutate(Name = ifelse(is.na(Name), CPD_NAME, Name))

cdrp.profiles <- cdrp.profiles %>%
  left_join(., moa2, by = "Metadata_broad_sample") %>%
  filter(!is.na(Name) & !is.na(MOA))

cdrp.profiles <- cdrp.profiles %>% 
  group_by(Metadata_broad_sample) %>%
  slice(1) %>%
  ungroup()

cdrp.profiles <- data.frame(cdrp.profiles)

rownames(cdrp.profiles) <- cdrp.profiles$Metadata_broad_sample

same.moa <- function(x, y) {
  if (is.na(x) || is.na(y) || x == "" || y == "") 
    return(FALSE)
  xs <- strsplit(x, ", ")[[1]]
  ys <- strsplit(y, ", ")[[1]]
  return(any(xs %in% ys) | any(ys %in% xs))
}

same.moa <- Vectorize(same.moa)

cr.melt <- cdrp.profiles %>%
  select(one_of(feats)) %>%
  t %>%
  cor %>% 
  reshape2::melt() %>%
  left_join(., cdrp.profiles %>% select(-one_of(feats)), by = c("Var1" = "Metadata_broad_sample")) %>%
  left_join(., cdrp.profiles %>% select(-one_of(feats)), by = c("Var2" = "Metadata_broad_sample")) %>%
  mutate(same.moa = same.moa(MOA.x, MOA.y))

same.val <- cr.melt %>% 
  filter(Var1 != Var2 & same.moa == 1) %>% 
  select(value) %>%
  as.matrix() %>%
  as.vector()

diff.val <- cr.melt %>% 
  filter(Var1 != Var2 & same.moa == 0) %>% 
  select(value) %>%
  as.matrix() %>%
  as.vector()

df1 <- data.frame(same.moa = "same MoA", correlation = (same.val)) 
df2 <- data.frame(same.moa = "different MoAs", correlation = (diff.val))

library(ggplot2)

g <- ggplot2::ggplot(mapping = aes(x = correlation, y=..density..)) +
  geom_histogram(data = df1, aes(fill = "red"), bins = 40, alpha=0.5) +
  geom_histogram(data = df2, aes(fill = "blue"), bins = 40, alpha=0.5) +
  xlab("Correlation") + 
  scale_fill_manual(name = "", values = c("red", "blue"), labels=c("red" = "same MoA", "blue" = "different MoA")) 
  
ggsave("same_vs_diff_moa_corr.png", g)
