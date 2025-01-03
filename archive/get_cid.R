all.digits <- function(str) {
  x <- tryCatch({as.numeric(str)}, warning = function(w) {NA})
  return(!is.na(x) & x == str)
}

cmpd.name.to.cid <- function(name) {
  v <- RCurl::getURL(sprintf("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/%s/cids/TXT", name))
  if (str_length(v) < 2) {
    return("")
  }
  a <- (str_sub(v, 1, (str_length(v) - 1)))
  b <- str_split(a, "\n")[[1]][1]
  if (!all.digits(b)) {
    return("")
  }
  return(b)
}

if (!file.exists("../results/manual/cmpd_cid.rds")) {
  cmpds <- Pf.cmpd$data$Image_Metadata_SOURCE_COMPOUND_NAME
  cmpd.cid <- data.frame(Image_Metadata_SOURCE_COMPOUND_NAME = c(), PubMedIDs = c())
  a <- progress::progress_bar$new(total = length(cmpds))

  for (cmpd in cmpds) {
    x <- -1
    while(x == -1) {
      x <- tryCatch({cmpd.name.to.cid(cmpd %>% str_to_lower())}, error = function(e) {-1})
    }
    cmpd.cid <- rbind(cmpd.cid, data.frame(Image_Metadata_SOURCE_COMPOUND_NAME = cmpd, PubMedIDs = x))
    a$tick()
  }

  saveRDS(cmpd.cid, "cmpd_cid.rds")
} else {
  cmpd.cid <- readRDS("../results/manual/cmpd_cid.rds")
}