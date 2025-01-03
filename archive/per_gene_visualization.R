rm(list = ls())

library(dplyr)
library(stringr)

load("../results/master/2017-12-06_9bc9b5ae_good/workspace.RData")

sign.chr <- Vectorize(function(x) ifelse(x > 0, "+", "-"))
hist.chr <- Vectorize(function(x) {
                        y <- unlist(str_split(x, ", "))
                        lvls <- unique(y)
                        hist.vec <- rep(0, length(lvls))
                        names(hist.vec) <- lvls

                        for (i in 1:length(y)) {
                          hist.vec[y[i]] <- hist.vec[y[i]] + 1
                        }

                        hist.vec <- sort(hist.vec, decreasing = T)
                        paste0(names(hist.vec), " (", unname(hist.vec), ")", collapse = ", ")
                      })

cmpd.annot <- gene.compound.cr.pr %>%
   group_by(Var2) %>%
   slice(1:k) %>%
   filter(valid >= 1) %>%
   ungroup() %>%
   arrange(-valid) %>%
   mutate(Var1 = paste0(Var1, " (", valid, ", ", sign.chr(value), ")")) %>%
   mutate(MOA = paste0(MOA, " (", valid, ")")) %>%
   group_by(Var2) %>%
   summarise(Name = paste0(unique(Var1), collapse = ", "),
             MOA = paste0(unique(MOA), collapse = ", "),
             Target = hist.chr(paste0(Target, collapse = ", ")))


data.frame(gene = .gene, p.value = p.vals) %>%
  filter(!gene %in% gene.blacklist) %>%
  mutate(p.value = p.adjust(p.value, "BH")) %>%
  arrange(p.value) %>%
  mutate(p.value = round(p.value, 3)) %>%
  left_join(., pthw, by = c("gene" = "Gene")) %>%
  rename(adjusted.p.value = p.value) %>%
  left_join(., cmpd.annot, by = c("gene" = "Var2")) %>%
  rename(Gene = gene, Gene.Pathway = Pathway, Compound.Name = Name) %>%
  htmlTable::htmlTable()

data.frame(gene = .gene, p.value = p.vals) %>%
  filter(!gene %in% gene.blacklist) %>%
  left_join(., Pf_org$data %>%
              select(Gene, Pathway) %>%
              unique, by = c("gene" = "Gene")) %>%
  mutate(p.value = round(p.adjust(p.value, "BH"), 4)) %>%
  arrange(p.value) %>%
  rename(Gene = gene) %>%
  select(Gene, Pathway, p.value) %>%
  htmlTable::htmlTable()
