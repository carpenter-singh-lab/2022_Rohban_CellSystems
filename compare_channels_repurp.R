rm(list = ls())
library(dplyr)
library(stringr)
library(ggplot2)

# DNA <- readRDS("../results/master/2017-06-09_de0bf46f/MOA_consistency__DNA.rds")
# RNA <- readRDS("../results/master/2017-06-09_de0bf46f/MOA_consistency__RNA.rds")
# Mito <- readRDS("../results/master/2017-06-09_de0bf46f/MOA_consistency__Mito.rds")
# ER <- readRDS("../results/master/2017-06-09_de0bf46f/MOA_consistency__ER.rds")
# AGP <- readRDS("../results/master/2017-06-09_de0bf46f/MOA_consistency__AGP.rds")

DNA <- readRDS("../results/master/2017-06-08_6cf001b8/MOA_consistency_DNA.rds")
RNA <- readRDS("../results/master/2017-06-08_6cf001b8/MOA_consistency_RNA.rds")
Mito <- readRDS("../results/master/2017-06-08_6cf001b8/MOA_consistency_Mito.rds")
ER <- readRDS("../results/master/2017-06-08_6cf001b8/MOA_consistency_ER.rds")
AGP <- readRDS("../results/master/2017-06-08_6cf001b8/MOA_consistency_AGP.rds")
AreaShape <- readRDS("../results/master/2017-06-12_149e5efa/MOA_consistency_AreaShape.rds")


lst <- list("DNA" = DNA, 
            "RNA" = RNA, 
            "Mito" = Mito, 
            "ER" = ER, 
            "AGP" = AGP,
            "AreaShape" = AreaShape)

chn <- names(lst)
agg <- NULL

for (ch in chn) {
  cl <- colnames(lst[[ch]])
  cl.new <- paste(cl, ch, sep = "_")
  cl.new[1] <- cl[1]
  colnames(lst[[ch]]) <- cl.new
  
  if (is.null(agg)) {
    agg <- lst[[ch]]
  } else {
    agg <- agg %>% left_join(., lst[[ch]], by = "MOA")
  }
}

agg <- agg %>% mutate(how.many.chnls = consistent_AGP + 
                        consistent_DNA + 
                        consistent_ER +
                        consistent_Mito +
                        consistent_RNA + 
                        consistent_AreaShape)

agg <- agg %>% dplyr::mutate(consistent_AGP = ifelse(is.na(consistent_AGP), 
                                                     FALSE,
                                                     consistent_AGP))

agg <- agg %>% dplyr::mutate(consistent_AreaShape = ifelse(is.na(consistent_AreaShape), 
                                                     FALSE,
                                                     consistent_AreaShape))

agg <- agg %>% dplyr::mutate(consistent_ER = ifelse(is.na(consistent_ER), 
                                                     FALSE,
                                                     consistent_ER))
agg <- agg %>% dplyr::mutate(consistent_Mito = ifelse(is.na(consistent_Mito), 
                                                     FALSE,
                                                     consistent_Mito))
agg <- agg %>% dplyr::mutate(consistent_DNA = ifelse(is.na(consistent_DNA), 
                                                     FALSE,
                                                     consistent_DNA))
agg <- agg %>% dplyr::mutate(consistent_RNA = ifelse(is.na(consistent_RNA), 
                                                     FALSE,
                                                     consistent_RNA))


agg1 <- agg %>% dplyr::group_by(MOA, how.many.chnls) %>%
  dplyr::summarise(which.chnls = paste(c("AGP", "DNA", "RNA", "Mito", "ER", "AreaShape")[c(consistent_AGP,
                                                             consistent_DNA, 
                                                             consistent_RNA, 
                                                             consistent_Mito,
                                                             consistent_ER,
                                                             consistent_AreaShape)],
                            collapse = "_"))  %>%
  filter(how.many.chnls > 0) %>%
  arrange(how.many.chnls) 

agg1 %>% htmlTable::htmlTable()

ggplot(agg, aes(x = how.many.chnls)) + 
  geom_histogram(aes(y=..count../sum(..count..))) +
  ylab("probability") + 
  xlab("Number of channels compounds in an MOA found to be correlated")

lvl <- agg1 %>% ungroup() %>% arrange(how.many.chnls) %>% dplyr::select(which.chnls, -MOA) %>% as.matrix() %>% as.vector() %>% unique 

ggplot(agg1 %>% mutate(which.chnls = factor(which.chnls, 
                                            levels = rev(lvl),
                                            ordered = T)), aes(x = which.chnls)) + 
  geom_histogram(aes(y=..count../sum(..count..)), stat="count") +
  ylab("probability") + coord_flip()

agg1 %>% 
  left_join(., agg, by = "MOA") %>%
  htmlTable::htmlTable()

agg2 <- agg1 %>% 
  left_join(., agg, by = "MOA") %>% 
  mutate(strongest = ifelse(which.chnls == "ER", sprintf("%f - %f", sig.strn_ER, thresh_ER),
         ifelse(which.chnls == "RNA", sprintf("%f - %f", sig.strn_RNA, thresh_RNA),
         ifelse(which.chnls == "AreaShape", sprintf("%f - %f", sig.strn_AreaShape, thresh_AreaShape),
         ifelse(which.chnls == "DNA", sprintf("%f - %f", sig.strn_DNA, thresh_DNA), 
         ifelse(which.chnls == "Mito", sprintf("%f - %f", sig.strn_Mito, thresh_Mito),
         ifelse(which.chnls == "AGP", sprintf("%f - %f", sig.strn_AGP, thresh_AGP),
                "")))))))

agg2 <- agg2 %>% 
  dplyr::select(MOA, how.many.chnls.x, which.chnls, n.members_DNA, strongest) %>%
  dplyr::rename(how.many.chnls = how.many.chnls.x, 
         n.members = n.members_DNA)

data.frame(frequency = apply(agg[which(agg$how.many.chnls <= 4),c("consistent_DNA",
             "consistent_RNA",
             "consistent_Mito",
             "consistent_ER",
             "consistent_AGP", 
             "consistent_AreaShape")], 2,
      sum)) %>% 
  tibble::rownames_to_column(var = "Channel") %>% 
  mutate(Channel = str_replace_all(Channel, "consistent_", "")) %>%
  arrange(-frequency) %>%
  mutate(Channel = factor(Channel, levels = Channel)) %>%
  ggplot(., aes(x = Channel, y = frequency)) + 
  geom_bar(stat = "identity", width = 0.2) +
  ylab("Number of MOAs with a strong signature \n in less than 5 channels") +
  xlab("Contributing channels")

data.frame(frequency = apply(agg[which(agg$how.many.chnls == 1),c("consistent_DNA",
                                                                  "consistent_RNA",
                                                                  "consistent_Mito",
                                                                  "consistent_ER",
                                                                  "consistent_AGP",
                                                                  "consistent_AreaShape")], 2,
                             sum)) %>% 
  tibble::rownames_to_column(var = "Channel") %>% 
  mutate(Channel = str_replace_all(Channel, "consistent_", "")) %>%
  arrange(-frequency) %>%
  mutate(Channel = factor(Channel, levels = Channel)) %>%
  ggplot(., aes(x = Channel, y = frequency)) + 
  geom_bar(stat = "identity", width = 0.2) +
  ylab("Number of MOAs with a strong signature \n in just one channel") +
  xlab("Contributing channels")
