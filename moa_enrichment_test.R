rm(list = ls())
library(dplyr)
library(stringr)
library(plyr)
library(reshape2)
library(htmlTable)
library(stringi)
library(foreach)
library(doMC)

doMC::registerDoMC(cores = 3)
set.seed(42)

source("rep.corr.func.R")
use.repurp.annots <- F
permute.moas <- F
profile.thresholding <- F
profile.thr <- 3
seed.moa <- -8     ## ignore, if permute.moas is False

base.dir <- "2016-07-20_90b7bb86"
use.feat.selected.in.TA.ORF.in.common <- F
corr.type <- "pearson"

x <- readRDS("../input/repurp/2016_04_01_a549_48hr_batch1_normalized.rds")
moa.list <- x[,c("Metadata_pert_iname", "Metadata_moa")]
m.list <- moa.list$Metadata_moa %>% unique %>% as.matrix() %>% as.vector() %>% as.vector()
m.list <- setdiff(m.list, NA)
set.seed(seed.moa)
pr <- permute::shuffle(length(m.list))
m.list <- data.frame(moa = m.list)
rownames(m.list) <- m.list$moa[pr]
m.list.p <- lapply(x$Metadata_moa, function(x) m.list[(x %>% as.character()), "moa"]) %>% unlist
if (permute.moas) {
  x$Metadata_moa <- m.list.p  
}

x <- cbind(x, data.frame(Metadata_Treatment = paste(x$Metadata_pert_id, x$Metadata_mg_per_ml, sep = "@")))
feats <- colnames(x)
feats <- feats[which(!str_detect(feats, "Metadata"))]
metadata <- colnames(x)
metadata <- metadata[which(str_detect(metadata, "Metadata"))]
if (profile.thresholding) {
  z <- x[,feats]
  z[z > profile.thr] <- profile.thr
  z[z < -profile.thr] <- -profile.thr
  x[,feats] <- z
}
thr <- non.rep.cor(list(data = x, feat_col = feats, factor_col = metadata), "Metadata_Treatment", feats)
u <- rep.cor(list(data = x, feat_col = feats, factor_col = metadata), "Metadata_Treatment", feats)
strong.trt <- u$Metadata_Treatment[which(u$cr > thr)]
sprintf("Hit ratio (compound-concentrations) : %f%%", round(length(strong.trt)/NROW(u) * 100))
strong.cmpd <- lapply(strong.trt, function(x) str_split(x, "@")[[1]][1]) %>% unlist %>% unique
all.cmpd <- lapply(u$Metadata_Treatment, function(x) str_split(x, "@")[[1]][1]) %>% unlist %>% unique
sprintf("Hit ratio (compounds) : %f%%", round(length(strong.cmpd)/length(all.cmpd) * 100))
x.all <- x
x <- x %>% dplyr::filter(Metadata_Treatment %in% strong.trt)
x.collapsed <- x %>% dplyr::group_by(Metadata_pert_iname, Metadata_pert_idose, Metadata_moa) %>% 
  dplyr::select(one_of(c(feats, "Metadata_pert_iname", "Metadata_pert_idose", "Metadata_moa"))) %>% dplyr::summarise_each(funs("mean")) %>% dplyr::ungroup() %>% dplyr::filter(!is.na(Metadata_pert_iname))

Pf.cmpd <- list(data = x.collapsed, feat_cols = feats, factor_cols = c("Metadata_pert_iname", "Metadata_pert_idose", "Metadata_moa"))
Pf.cmpd.all <- list(data = x.all, feat_cols = feats, factor_cols = c("Metadata_pert_iname", "Metadata_pert_idose", "Metadata_moa"))

strong.cmpd <- Pf.cmpd$data %>% dplyr::group_by(Metadata_pert_iname) %>% dplyr::summarise(cnt = n()) %>%
  dplyr::filter(cnt >= 4) %>% dplyr::select(Metadata_pert_iname) %>% unique

Pf.cmpd$data <- Pf.cmpd$data %>% dplyr::filter(Metadata_pert_iname %in% (strong.cmpd %>% as.matrix() %>% as.vector()))

distinct.moas <- lapply(Pf.cmpd$data$Metadata_moa, function(x) (str_split(x, "\\|")[[1]])) %>% unlist %>% 
  unique %>% setdiff(., NA)

agg.fn <- function(x) return(ifelse(quantile(x, 0.7) > -quantile(x, 0.3), quantile(x, 0.7), quantile(x, 0.30)))
#agg.fn <- mean

cons <- data.frame(MOA = c(), consistent = c(), sig.strn = c(), thresh = c(), n.members = c())
n.sample <- 100

cons <- foreach(moa = distinct.moas) %dopar% {
  Px <- Pf.cmpd$data %>% dplyr::filter(stri_detect_fixed(Metadata_moa, moa)) 
  cr <- cor(Px[,Pf.cmpd$feat_cols] %>% t, method = corr.type)
  colnames(cr) <- Px$Metadata_pert_iname
  rownames(cr) <- Px$Metadata_pert_iname
  
  cr1 <- cr %>% reshape2::melt() %>% 
    dplyr::group_by(Var1, Var2) %>% dplyr::summarise(agg.value = agg.fn(value)) %>%
    reshape2::acast(Var1 ~ Var2, value.var = "agg.value")
  
  smp <- c()
  
  i <- 1
  while (i <= n.sample) {
    cmpd.sm <- sample(Pf.cmpd$data$Metadata_pert_iname %>% unique, NROW(cr1))
    Px <- Pf.cmpd$data %>% dplyr::filter(Metadata_pert_iname %in% cmpd.sm) 
    if (length(unique(setdiff(Px$Metadata_moa, NA))) < length((setdiff(Px$Metadata_moa, NA)))) {
      next
    }
    i <- i + 1
    
    cr <- cor(Px[,Pf.cmpd$feat_cols] %>% t, method = corr.type)
    cr1.tmp <- cr %>% reshape2::melt() %>% 
      dplyr::group_by(Var1, Var2) %>% dplyr::summarise(agg.value = agg.fn(value)) %>%
      reshape2::acast(Var1 ~ Var2, value.var = "agg.value")
    smp <- c(smp, cr1.tmp %>% as.dist() %>% median(., na.rm = T))
  }
  thr <- quantile(smp, 0.95, na.rm = T)

  v <- cr1 %>% as.dist() %>% as.vector() %>% median
  data.frame(MOA = moa, consistent = (v >= thr), sig.strn = round(v, 2),
             thresh = round(thr, 2), n.members = NROW(cr1))
}

MOA.consistency <- do.call(rbind, cons)
MOA.consistency %>% dplyr::arrange(-(sig.strn - thresh)) %>% htmlTable()
sum(MOA.consistency$consistent, na.rm = T)/(length(which(MOA.consistency$n.members > 1)))

