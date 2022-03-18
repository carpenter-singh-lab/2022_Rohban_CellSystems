## target test

load("../input/TA/Initial_analysis_workspace_new.RData")
moa2 <- readr::read_csv("../input/CDP2/MOA_annot2.csv")

all.targ <- lapply(moa2$Target %>% 
                     unique, 
                   function(x) str_split(x, ", ")[[1]]) %>% 
  do.call(c, .) %>% 
  unique

path <- sprintf("../results/manual/ppis")

if (!dir.exists(path)) {
  system(sprintf("mkdir -p %s", path))
}

if (file.exists(sprintf("%s/all.rds", path))) {
  ppi.all <- readRDS(sprintf("%s/all.rds", path))
} else {
  ppi.all <- tryCatch(get.all.interacting.proteins(all.targ), error = function(e) {})
  saveRDS(ppi.all, sprintf("%s/all.rds", path))
}

path <- sprintf("../results/manual/ppis")
ppi.all <- readRDS(sprintf("%s/all.rds", path))

genes <- Pf.strong$Gene[str_detect(Pf.strong$AlleleDesc, "WT")]

ppi.all <- ppi.all %>%
  filter(Protein.1 %in% genes | Protein.2 %in% genes)

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

valid.pair <- function(gene.name, compound.targ) {
  target <- compound.targ
  
  if (length(target) == 0) {
    return(F)
  }
  
  if (is.null(target)) {
    return(F)
  }
  
  if (is.na(target)) {
    return(F)
  }
  
  if(target == "") {
    return(F)
  }
  
  targets <- str_split(target, ", ")[[1]]
  
  if (gene.name %in% targets) {
    return(T)
  }
  
  n <- ppi.all %>% 
    dplyr::filter(Protein.1 == gene.name & Protein.2 %in% targets |
                    Protein.2 == gene.name & Protein.1 %in% targets) %>%
    NROW
  
  return(n)
}

library(ontologyIndex)
data(go)

library(ontologySimilarity)
data(gene_GO_terms)
data(GO_IC)

valid.pair.alt <- function(gene.name, compound.targ) {
  target <- compound.targ
  
  if (length(target) == 0) {
    return(0)
  }
  
  if (is.null(target)) {
    return(0)
  }
  
  if (is.na(target)) {
    return(0)
  }
  
  if(target == "") {
    return(0)
  }
  
  targets <- str_split(target, ", ")[[1]]

  beach <- gene_GO_terms[c(gene.name, targets)]
  
  beach <- beach[which(!unlist(lapply(beach, is.null)))]
  
  if (length(beach) == 1) {
    return(0)
  }
  
  cc <- go$id[go$name == "biological_process"]
  #print(beach)
  beach_cc <- lapply(beach, function(x) intersection_with_descendants(go, roots=cc, x)) 
  data.frame(check.names=FALSE, `#terms`=sapply(beach, length), `#CC terms`=sapply(beach_cc, length))
  
  sim_matrix_cc <- get_sim_grid(
    ontology=go, 
    information_content=GO_IC,
    term_sets=beach_cc)
  
  tg <- colnames(sim_matrix_cc)
  n <- max(sim_matrix_cc[gene.name, tg[2:length(tg)]])
  n <- as.vector(as.matrix(unname(n)))
  
  return(n)
}
