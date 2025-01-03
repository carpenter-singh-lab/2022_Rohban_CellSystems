library(dplyr)
library(stringr)

source("read_dataset.R")
source("rep.corr.func.R")

Pf <- read.dataset("Repurposing", dose.closest = NULL)

cmpds <- readr::read_csv("../input/dopamine_drugs.csv", col_names = F)

v <- non.rep.cor(Pf, "Metadata_Treatment", Pf$feat_cols)
u <- rep.cor(Pf, "Metadata_Treatment", Pf$feat_cols)

df <- Pf$data %>%
  filter(str_sub(Metadata_broad_sample, 1, 13) %in% str_sub(cmpds$X1, 1, 13))

u <- u %>%
  filter(Metadata_Treatment %in% df$Metadata_Treatment) %>%
  arrange(-cr) %>%
  left_join(., Pf$data[, c("Metadata_pert_iname",
                           "Metadata_mmoles_per_liter",
                           "Metadata_Treatment")] %>% unique,
            by = "Metadata_Treatment") %>%
  arrange(-cr) %>%
  ungroup() %>%
  select(Metadata_Treatment, Metadata_pert_iname, Metadata_mmoles_per_liter, cr)

cmpd.strong <- u %>%
  filter(cr > v) %>%
  select(Metadata_Treatment) %>%
  as.matrix() %>%
  as.vector()

metadata <- c("Metadata_pert_iname",
              "Metadata_broad_sample",
              "Metadata_mmoles_per_liter",
              "Metadata_moa")

df <- df %>%
  filter(Metadata_Treatment %in% cmpd.strong) %>%
  select(one_of(c(Pf$feat_cols,
                  metadata))) %>%
  group_by_at(.vars = metadata) %>%
  summarise_all(.funs = "mean") %>%
  ungroup()

cr <- df %>%
  select(one_of(Pf$feat_cols)) %>%
  t %>%
  cor

df <- df %>%
  mutate(Treatment = paste0(Metadata_pert_iname,
                            " @ ",
                            Metadata_mmoles_per_liter))

colnames(cr) <- df$Treatment
rownames(cr) <- df$Treatment

cr %>%
  reshape2::melt() %>%
  left_join(., df[, c(metadata, "Treatment")], by = c("Var1" = "Treatment")) %>%
  left_join(., df[, c(metadata, "Treatment")], by = c("Var2" = "Treatment")) %>%
  arrange(-value) %>%
  select(Metadata_pert_iname.x,
         Metadata_pert_iname.y,
         Metadata_mmoles_per_liter.x,
         Metadata_mmoles_per_liter.y,
         Metadata_moa.x,
         Metadata_moa.y,
         value) %>%
  group_by(Metadata_pert_iname.x, Metadata_pert_iname.y) %>%
  summarise(value = max(value)) %>%
  ungroup() %>%
  reshape2::acast(Metadata_pert_iname.x ~ Metadata_pert_iname.y) %>%
  corrplot::corrplot(.,
                     order = "hclust",
                     hclust.method = "average")
