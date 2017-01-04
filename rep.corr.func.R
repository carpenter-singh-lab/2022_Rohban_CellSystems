library(dplyr)
library(reshape2)
library(permute)

rep.cor <- function(Pf, grp.var, feat.var) {
  x <- Pf$data %>% dplyr::select(one_of(c(grp.var, feat.var))) %>%
    dplyr::group_by_(grp.var) %>% do(data.frame(cr = median(as.dist(cor(t(.[,feat.var]))), na.rm = T)))
  return(x)
}

non.rep.cor <- function(Pf, grp.var, feat.var) {
  pr <- permute::shuffle(NROW(Pf$data))
  Pf$data[,grp.var] <- Pf$data[pr,grp.var]
  u <- rep.cor(Pf, grp.var, feat.var)
  return(quantile(u$cr, 0.95, na.rm = T))
}