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