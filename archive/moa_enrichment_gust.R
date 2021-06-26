rm(list = ls())
library(dplyr)
library(stringr)
library(plyr)
library(reshape2)
library(htmlTable)
library(stringi)
library(foreach)
library(doMC)

doMC::registerDoMC(cores = 2)

source("rep.corr.func.R")

use.repurp.annots <- F
permute.moas <- F
seed.moa <- -8     ## ignore, if permute.moas is False

base.dir <- "2016-07-20_90b7bb86"
use.feat.selected.in.TA.ORF.in.common <- F
corr.type <- "pearson"

base.dir <- "2016-07-20_90b7bb86"
cmpd.name <- readRDS(sprintf("../input/Gustafsdottir/%s/cmpds.names.rds", base.dir))

load("../input/Gustafsdottir/Initial_analysis.RData")
u <- rep.cor(Pf.full.plate.norm, "Image_Metadata_BROAD_ID", Pf.full.plate.norm$feat_cols)
thr <- non.rep.cor(Pf.full.plate.norm, "Image_Metadata_BROAD_ID", Pf.full.plate.norm$feat_cols)
strong.cmpd <- u$Image_Metadata_BROAD_ID[which(u$cr > thr)]

Pf.final <- Pf.full.plate.norm$data %>% dplyr::filter(Image_Metadata_BROAD_ID %in% strong.cmpd) %>% as.data.frame() 
Pf.final.all <- Pf.full.plate.norm$data %>% as.data.frame() 

all.cols <- colnames(Pf.final)
fact.col <- all.cols[which(str_detect(all.cols, "Metadata") | all.cols %in% c("Plate", "Well"))]
feat.col <- setdiff(all.cols, fact.col)
Pf.gust <- list(data = Pf.final, feat_cols = feat.col, factor_cols = fact.col)
Pf.gust.all <- list(data = Pf.final.all, feat_cols = feat.col, factor_cols = fact.col)

Pf.gust$data <- Pf.gust$data %>% dplyr::select(one_of(c("Image_Metadata_BROAD_ID", Pf.gust$feat_cols))) %>% dplyr::group_by(Image_Metadata_BROAD_ID) %>% dplyr::summarise_each(funs("mean"))

Pf.gust$data <- plyr::join(Pf.gust$data, cmpd.name, by = "Image_Metadata_BROAD_ID")
Pf.gust$factor_cols <- c(Pf.gust$factor_cols, "Image_Metadata_SOURCE_COMPOUND_NAME")
Pf.cmpd <- Pf.gust
Pf.cmpd.all <- Pf.gust.all
Pf.cmpd.all$data <- Pf.cmpd.all$data[,colnames(Pf.cmpd$data)]
Pf.cmpd.all$feat_cols <- Pf.cmpd$feat_cols
Pf.cmpd.all$factor_cols <- Pf.cmpd$factor_cols
Pf.cmpd$data$Image_Metadata_SOURCE_COMPOUND_NAME <- lapply(Pf.cmpd$data$Image_Metadata_SOURCE_COMPOUND_NAME, 
                                                           function(x) str_to_lower(x)) %>% unlist
moas <- read.csv("../input/cmpd_moa.csv", header = T)
moas$Name <- lapply(moas$Name, function(x) str_to_lower(x)) %>% unlist
colnames(moas)[1] <- "Image_Metadata_SOURCE_COMPOUND_NAME"
Pf.cmpd$data <- plyr::join(Pf.cmpd$data, moas[,c("Image_Metadata_SOURCE_COMPOUND_NAME", "MOA")], by = "Image_Metadata_SOURCE_COMPOUND_NAME")
Pf.cmpd$factor_cols <- c(Pf.cmpd$factor_cols, "MOA")
distinct.moas <- unique(Pf.cmpd$data$MOA)

MOA.consistency <- data.frame(MOA = c(), consistent = c(), sig.strn = c(), thresh = c(), n.members = c())
n.sample <- 20

cons <- foreach(moa = distinct.moas) %dopar% {
  Px <- Pf.cmpd$data %>% dplyr::filter(MOA == moa) 
  cr <- cor(Px[,Pf.cmpd$feat_cols] %>% t, method = corr.type)
  colnames(cr) <- Px$Image_Metadata_SOURCE_COMPOUND_NAME
  rownames(cr) <- Px$Image_Metadata_SOURCE_COMPOUND_NAME
  
  cr1 <- cr 
  
  smp <- c()
  
  i <- 1
  while (i <= n.sample) {
    cmpd.sm <- sample(Pf.cmpd$data$Image_Metadata_SOURCE_COMPOUND_NAME %>% unique, NROW(cr1))
    Px <- Pf.cmpd$data %>% dplyr::filter(Image_Metadata_SOURCE_COMPOUND_NAME %in% cmpd.sm) 
    if (length(unique(setdiff(Px$MOA, NA))) < length((setdiff(Px$MOA, NA)))) {
      next
    }
    i <- i + 1
    cr <- cor(Px[,Pf.cmpd$feat_cols] %>% t, method = corr.type)
    cr1.tmp <- cr 
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

