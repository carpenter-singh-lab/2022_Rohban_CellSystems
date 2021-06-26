library(dplyr)
library(stringr)

cdrp1 <- readRDS("../results/master/2017-04-20_425d653/Pf_bio_new.rds")
cdrp2 <- readRDS("../results/master/2017-04-20_425d653/Pf_DOS_new.rds")
presel <- c("leflunomide", "spironolactone", "nitrofural")
cdrp <- rbind(cdrp1, cdrp2)

feats <- colnames(cdrp)
feats <- feats[which(!str_detect(feats, "Metadata_"))]

all.cons <- readRDS("../results/master/2017-09-16_e1d600ab/gene_compound_all.rds")
cmpds <- all.cons %>% 
  filter(gene == "SMAD3_WT.1" & specificity > 0.1 & !CPD_NAME %in% presel) %>%
  select(Metadata_broad_sample) %>%
  as.matrix() %>%
  as.vector()

cdrp.sub <- cdrp %>% 
  filter(Metadata_broad_sample %in% cmpds) %>%
  group_by(Metadata_broad_sample) %>%
  summarise_at(.vars = feats, .funs = mean)
  
cr <- cor(cdrp.sub %>% select(-Metadata_broad_sample) %>% t)
rownames(cr) <- cdrp.sub$Metadata_broad_sample
colnames(cr) <- cdrp.sub$Metadata_broad_sample

hcl <- hclust(as.dist(1 - cr), method = "complete")
ct <- cutree(hcl, k = 7)

clst <- ct %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column("Metadata_broad_sample") %>%
  rename(!!"cluster" := !!".")

grp1 <- all.cons %>%
  filter(gene == "SMAD3_WT.1" & specificity > 0.1 & !CPD_NAME %in% presel) %>%
  left_join(., clst) %>%
  arrange(-specificity) %>%
  group_by(cluster) %>%
  slice(1) %>%
  ungroup() %>% 
  select(-cluster)

grp2 <- all.cons %>%
  filter(gene == "SMAD3_WT.1" & CPD_NAME %in% presel) 

cmpd <- rbind(grp1, grp2)

readr::write_csv(cmpd, "SMAD_cmpd.csv")