library(dplyr)
library(stringr)

a <- read.csv("../input/profiles.csv")
query.subpops <- function(enr.subpops, denr.subpops) {
  cmpd <- a$compound %>% unique
  for (x in enr.subpops) {
    y <- a %>% dplyr::filter(str_detect(subpop.enr, sprintf("_%s_", x))) %>%
      dplyr::select(compound) %>% as.matrix() %>% as.vector()
    cmpd <- intersect(cmpd, y)
  }
  for (x in denr.subpops) {
    y <- a %>% dplyr::filter(str_detect(subpop.denr, sprintf("_%s_", x))) %>%
      dplyr::select(compound) %>% as.matrix() %>% as.vector()
    cmpd <- intersect(cmpd, y)
  }
  return(cmpd)
}

query.subpops(enr.subpops = c("slightly elongated", "large and symmetric org. high DNA content"), denr.subpops = 
                c("asymmetric org. high DNA content", 
                  "small and asymmetric org. high DNA content"))

x <- readRDS("../input/repurp/2016_04_01_a549_48hr_batch1_normalized.rds")
x <- x %>% dplyr::filter(Metadata_moa == "glucocorticoid receptor agonist") %>% dplyr::select(Metadata_pert_iname) %>% unique
a %>% dplyr::filter(compound %in% (x %>% as.matrix() %>% as.vector())) %>% View 
