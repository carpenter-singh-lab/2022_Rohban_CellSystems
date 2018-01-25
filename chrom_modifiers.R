rm(list = ls())
library(dplyr)
library(stringr)
source("read_dataset.R")
source("rep.corr.func.R")

get.brd <- function(x) {
  paste(str_split(x, "-")[[1]][1:2], collapse = "-")
}

get.brd <- Vectorize(get.brd)

a <- read.csv("../input/results-5.txt", sep = "\t")
brds <- unique(a$Broad.ID)
brds.q <- unlist(lapply(brds, function(x) paste(str_split(x, "-")[[1]][1:2], collapse = "-")))

Pf <- read.dataset("Repurposing", dose.closest = 10)
Pf$data <- Pf$data %>% mutate(Metadata_broad_id = get.brd(Metadata_broad_sample))

u <- rep.cor(Pf, "Metadata_broad_id", Pf$feat_cols)
v <- non.rep.cor(Pf, "Metadata_broad_id", Pf$feat_cols, quant = 0.95)

Px <- Pf$data %>%
  group_by(Metadata_Treatment) %>% 
  summarise_at(.vars = Pf$feat_cols, .funs = "mean")

Px <- Px %>%
  left_join(Pf$data %>% 
              filter(Metadata_broad_sample != "DMSO") %>%
              select(Metadata_Treatment, Metadata_moa, Metadata_broad_sample, Metadata_pert_iname, Metadata_broad_id) %>% 
              unique, 
            by = "Metadata_Treatment")

cr <- cor(t(Px[, Pf$feat_cols]))
rownames(cr) <- Px$Metadata_broad_id
colnames(cr) <- Px$Metadata_broad_id

brd.interest <- intersect(brds.q, Px$Metadata_broad_id %>% unique)
cr <- cr[brd.interest, ]

cr.melt <- cr %>% reshape2::melt()

cr.melt <- cr.melt %>%
  left_join(Px %>% select(starts_with("Metadata")), by = c("Var1" = "Metadata_broad_id")) %>%
  left_join(Px %>% select(starts_with("Metadata")), by = c("Var2" = "Metadata_broad_id")) 

cr.melt <- cr.melt %>%
  filter(!Var2 %in% brds.q)

cr.melt %>%
  filter(value > 0.7) %>%
  group_by(Metadata_moa.y) %>%
  tally() %>%
  arrange(-n) %>% 
  slice(1:5) %>%
  knitr::kable()

Px.sub <- Px %>% filter(Metadata_broad_id %in% brd.interest) 
rownames(Px.sub) <- Px.sub$Metadata_pert_iname

Px.sub %>% 
  select(starts_with("Nuclei_Texture_")) %>% 
  select(contains("_DNA_")) %>% 
  View
