## process all the matches to create a consice list of genes to prioritize them

rm(list = ls())
library(dplyr)
library(stringr)

res.all <- readRDS("../results/master/2017-05-17_15389033/gene_compound_all.rds") ## bioactives
#res.all <- readRDS("../results/master/2017-05-17_e971110f/gene_compound_all.rds") ## all

split <- function(x, i) {
  str_split(x, "_")[[1]][i]  
}

split <- Vectorize(split)

get.interacting.proteins <- function(protein.name) {
  dr <- "../results/manual/ppis"
  if (!dir.exists(dr)) {
    system(sprintf("mkdir -p %s", dr))
  }
  f.name <- sprintf("%s/%s.rds", dr, protein.name)
  if (file.exists(f.name)) {
    return(readRDS(f.name))
  }
  
  tbl <- read.table(sprintf("http://webservice.thebiogrid.org/interactions/?searchNames=true&geneList=%s&taxId=9606&includeHeader=true&accesskey=ca950b072394ce1897811022f7757222", protein.name), sep="\t", header = FALSE, fill = T)
  tbl <- tbl[, c(8, 9, 12, 13, 14)] 
  colnames(tbl) <- c("Protein.1", "Protein.2", "Method", "Type", "Evidence")
  saveRDS(tbl, f.name)
  return(tbl)
}

get.all.interacting.proteins <- function(protein.name) {
  res <- c()
  for (p in protein.name) {
    ppi <- tryCatch(get.interacting.proteins(p), error = function(e) {})
    res <- rbind(res, ppi)
  }
  return(res)
}

interacts_with <- function(gene, targets) {
  if (is.na(targets)) {
    return(FALSE)
  }
  
  tg <- str_split(targets, ", ")[[1]]
  p1 <- get.all.interacting.proteins(tg)  
  p2 <- get.all.interacting.proteins(gene)
  p1 <- c(p1$Protein.1 %>% as.character(), p1$Protein.2 %>% as.character()) %>% unique
  p2 <- c(p2$Protein.1 %>% as.character(), p2$Protein.2 %>% as.character()) %>% unique
  
  return((gene %in% p1) | (any(tg %in% p2)))
}

interacts_with <- Vectorize(interacts_with)

res <- res.all %>% 
  dplyr::mutate(is.wt = str_detect(gene, "_WT")) %>% 
  dplyr::mutate(gene.name = split(gene, 1)) %>%
  dplyr::mutate(allele.name = split(gene, 2)) %>%
  dplyr::mutate(allele.name = ifelse(is.wt, "WT", allele.name)) %>%
  dplyr::mutate(validated = interacts_with(gene.name, Target))

res.sm <- res %>%
  dplyr::group_by(gene.name, allele.name) %>%
  dplyr::summarise(no.high.corr.in.l1000 = sum(l1000.score > 0.97, na.rm = T), 
                   avg.cell.count = mean(Mean_Cells_Count), 
                   avg.validated = sum(validated, na.rm = T),
                   max.CP.mimics.score = max(Corr.),
                   max.CP.antimimics.score = min(Corr.)) 

res.sm %>% dplyr::mutate(avg.cell.count = round(avg.cell.count, 2)) %>% htmlTable::htmlTable()

res %>% write.csv(., "all_cons.csv", row.names = F)