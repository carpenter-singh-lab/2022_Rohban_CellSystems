library(dplyr)
library(reshape2)
library(permute)

rep.cor <- function(Pf, grp.var, feat.var, aux.var = NULL) {
  x <- Pf$data %>% dplyr::select(one_of(c(grp.var, feat.var))) %>%
    dplyr::group_by_(grp.var) %>% do(data.frame(cr = median(as.dist(cor(t(.[,feat.var]))), na.rm = T)))
  y <- Pf$data %>% dplyr::select(one_of(c(grp.var, feat.var))) %>%
    #dplyr::group_by_(grp.var) %>% do(data.frame(strn = sum(abs(.[,feat.var]) > 1)/(NROW(.[,feat.var]) * 
    #                                                                                     length(feat.var))))
    dplyr::group_by_(grp.var) %>% do(data.frame(strn = 1))
    
  z <- left_join(x, y, by = grp.var)
  if (!is.null(aux.var)) {
    w <- Pf$data[,c(grp.var, aux.var)] %>% 
      unique() %>%
      left_join(., z, by = grp.var)
  } else {
    w <- z
  }
  w %<>% dplyr::mutate(mac = (cr * strn)^0.5)
  return(w)
}

non.rep.cor.rob <- function(Pf, grp.var, feat.var, quant = 0.95) {
  us <- c()
  for (i in 1:10) {
    pr <- permute::shuffle(NROW(Pf$data))
    Pf$data[,grp.var] <- Pf$data[pr,grp.var]
    u <- rep.cor(Pf, grp.var, feat.var)
    us <- rbind(us, u)
  }
  return(quantile(us$cr, quant, na.rm = T))
}

non.rep.cor <- function(Pf, grp.var, feat.var, quant = 0.95) {
  pr <- permute::shuffle(NROW(Pf$data))
  Pf$data[,grp.var] <- Pf$data[pr,grp.var]
  u <- rep.cor(Pf, grp.var, feat.var)
  return(quantile(u$cr, quant, na.rm = T))
}