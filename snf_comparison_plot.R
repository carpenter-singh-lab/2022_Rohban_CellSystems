library(ggplot2)
library(dplyr)

a <- readRDS("../results/master/2017-06-07_660cbe38/data_eval_prof_ord.rds")
b <- readRDS("../results/master/2017-06-07_660cbe38/data_eval_prof_snf.rds")

a <- rbind(a %>% mutate(type = "normal profile"), 
           b %>% mutate(type = "SNF"))

g <- ggplot(a, aes(x = thr, y = enr.ratio, 
                   color = type)) + geom_line() + 
  xlab("Corr. threshold quantile") +
  ylab("Enrichment ratio")

ggsave("SNF_vs_profile.png", g, width = 8, height = 5)
