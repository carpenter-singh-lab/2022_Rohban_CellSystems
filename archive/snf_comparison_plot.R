library(ggplot2)
library(dplyr)

a <- readRDS("../results/master/2017-06-07_660cbe38/data_eval_prof_ord.rds")
b <- readRDS("../results/master/2017-06-07_660cbe38/data_eval_prof_snf.rds")
c <- readRDS("../results/master/2017-06-07_ee0ef21a/data_eval_prof_snf.rds")
d <- readRDS("../results/master/2017-06-07_0b23130c/data_eval_prof.rds")

a <- rbind(a %>% mutate(type = "normal profile"),
           b %>% mutate(type = "SNF (t = 10)"),
           c %>% mutate(type = "SNF (t = 20)"),
           d %>% mutate(type = "SNF (t = 15)"))

g <- ggplot(a, aes(x = thr, y = enr.ratio,
                   color = type)) + geom_line() +
  xlab("Corr. threshold quantile") +
  ylab("Enrichment ratio")

ggsave("SNF_vs_profile.png", g, width = 8, height = 5)
