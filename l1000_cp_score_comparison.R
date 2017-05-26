## plotting L1000 vs CP scores for gene-compounds
rm(list = ls())
library(dplyr)
library(ggplot2)

cr.melt.cp <- readRDS("../results/master/2017-05-02_9796dc34/cr_melt_cp.rds")
cr.melt.l1000 <- readRDS("../results/master/2017-05-02_aaa7dd14/cr_melt_l1000.rds")

cr.jn <- cr.melt.cp %>% 
  dplyr::left_join(., cr.melt.l1000, by = c("Var1", "Var2"))

ggplot(cr.jn, aes(x = value.x, y = value.y)) + geom_point()