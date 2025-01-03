rm(list = ls())
library(dplyr)

brd.samples <- read.csv("../followups/hras/compounds/samples.txt", sep = "\t")

kras.wt <- read.csv("../followups/hras/compounds/KRAS_WT.csv")
kras.wt <- kras.wt %>% mutate(match.to = "KRAS_WT")

hras.g12v <- read.csv("../followups/hras/compounds/HRAS_G12V.csv")
hras.g12v <- hras.g12v %>% mutate(match.to = "HRAS_G12V")

kras.g12v <- read.csv("../followups/hras/compounds/KRAS_G12V.csv")
kras.g12v <- kras.g12v %>% mutate(match.to = "KRAS_G12V")
kras.g12v <- cbind(kras.g12v, data.frame(Comments = rep(NA, NROW(kras.g12v))))

wt.mut.diff <- read.csv("../followups/hras/compounds/KRAS_diff_MUT_WT.csv")
wt.mut.diff <- wt.mut.diff %>% mutate(match.to = "WT Diff Mutant")

all.brds <- read.csv("../followups/hras/compounds/Compound_List.csv")

x <- rbind(kras.wt, hras.g12v, kras.g12v, wt.mut.diff)

x <- x %>%
  group_by(Broad.ID, match.to) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(Comments = ifelse(is.na(Comments), "", as.character(Comments)))

y <- x %>%
  group_by(Broad.ID) %>%
  summarise(Corr..Sign = paste0(Corr..Sign, collapse = "; "),
            Compound.Type = ifelse(any(Compound.Type == "Bioactive"), "Bioactive", "DOS"),
            match.to = paste0(match.to, collapse = "; "),
            Compound.Name = Compound.Name[1],
            Comments = paste0(Comments, collapse = "")) %>%
  ungroup() %>%
  mutate(Broad.ID = factor(Broad.ID, levels = all.brds$Broad.ID)) %>%
  arrange(Broad.ID) %>%
  right_join(., brd.samples, by = "Broad.ID")

readr::write_csv(y, "../followups/hras/all_compounds_annot.csv")