rm(list = ls())
library(dplyr)
library(stringr)
library(htmlTable)
library(reshape2)

cdp2 <- F
gust <- F

if (cdp2) {
  Pf.cdp2 <- readRDS("../input/CDP2/Pf.rds")
  Pf.cdp2$data <- Pf.cdp2$data %>% dplyr::filter(str_detect(Image_Metadata_CPD_Plate_Map_Name, "BIOA"))
  med <- apply(Pf.cdp2$data[,Pf.cdp2$feat_cols], 2, median)
  md <- apply(Pf.cdp2$data[,Pf.cdp2$feat_cols], 2, mad)
  Pf.cdp2$data[,Pf.cdp2$feat_cols] <- scale(Pf.cdp2$data[,Pf.cdp2$feat_cols], center = med, scale = md)
  
  all.cmpd <- Pf.cdp2$data$Image_Metadata_BroadID %>% unique
  rep.cor.list <- data.frame(trt = c(), rep.cor = c())
  nulls <- c()
  
  for (cmpd in all.cmpd) {
    x <- Pf.cdp2$data %>% dplyr::filter(Image_Metadata_BroadID == cmpd)
    cr.med <- cor(x[,Pf.cdp2$feat_cols] %>% t) %>% as.dist %>% median
    cmpd.samp <- sample(all.cmpd, 2)
    y <- Pf.cdp2$data %>% dplyr::filter(Image_Metadata_BroadID %in% cmpd.samp)
    cr.med.null <- cor(y[,] %>% dplyr::group_by(Image_Metadata_BroadID) %>% slice(1) %>% dplyr::ungroup() %>% dplyr::select(one_of(Pf.cdp2$feat_cols)) %>% t) %>% as.dist
    nulls <- c(nulls, cr.med.null)  
    rep.cor.list <- rbind(rep.cor.list, data.frame(trt = cmpd, rep.cor = cr.med))
  }
  thr <- quantile(nulls, 0.95, na.rm = T)
  #thr <- -1
  trt <- rep.cor.list$trt[which(rep.cor.list$rep.cor > thr)]
  Pf.cdp2$data <- Pf.cdp2$data %>% dplyr::filter(Image_Metadata_BroadID %in% trt)
  Pf.cdp2$data <- Pf.cdp2$data %>% dplyr::select(one_of(c("Image_Metadata_BroadID", Pf.cdp2$feat_cols))) %>%
    dplyr::group_by(Image_Metadata_BroadID) %>% dplyr::summarise_each(funs("mean"))
  i <- which(colnames(Pf.cdp2$data) == "Image_Metadata_BroadID")
  colnames(Pf.cdp2$data)[i] <- "Image_Metadata_BROAD_ID"
  Pf.cdp2$factor_cols <- c("Image_Metadata_BROAD_ID", setdiff(Pf.cdp2$factor_cols, "Image_Metadata_BroadID"))
  
  Pf.cmpd <- Pf.cdp2
} else if (gust) {
  load("../input/Gustafsdottir/Initial_analysis.RData")
  Pf.final <- x %>% dplyr::filter(Image_Metadata_BROAD_ID %in% strong.cmpd) %>% as.data.frame()
  rownames(Pf.final) <- (Pf.final[,1] %>% as.matrix() %>% as.vector())
  Pf.gust <- list(data = Pf.final, feat_cols = colnames(Pf.final)[2:NCOL(Pf.final)], factor_cols = colnames(Pf.final)[1])
  base.dir <- "2016-07-20_90b7bb86"
  cmpd.name <- readRDS(sprintf("../input/Gustafsdottir/%s/cmpds.names.rds", base.dir))
  Pf.gust$data <- plyr::join(Pf.gust$data, cmpd.name, by = "Image_Metadata_BROAD_ID")
  Pf.gust$factor_cols <- c(Pf.gust$factor_cols, "Image_Metadata_SOURCE_COMPOUND_NAME")
  Pf.cmpd <- Pf.gust
  med <- apply(Pf.cmpd$data[,Pf.cmpd$feat_cols], 2, median)
  md <- apply(Pf.cmpd$data[,Pf.cmpd$feat_cols], 2, mad)
  Pf.cmpd$data[,Pf.cmpd$feat_cols] <- scale(Pf.cmpd$data[,Pf.cmpd$feat_cols], med, md)
} else {
  Pf.cmpd <- readRDS("../input/CDP2/subset/Pf.rds")
  strong.cmpd <- readRDS("../input/CDP2/subset/hits_mean.profile.rds")
  Pf.cmpd$data <- Pf.cmpd$data %>% dplyr::filter(Image_Metadata_BROAD_ID %in% strong.cmpd)
  Pf.cmpd$data <- Pf.cmpd$data %>% dplyr::select(one_of(c("Image_Metadata_BROAD_ID", Pf.cmpd$feat_cols))) %>%
    dplyr::group_by(Image_Metadata_BROAD_ID) %>% dplyr::summarise_each(funs("mean"))
  base.dir <- "2016-07-20_90b7bb86"
  cmpd.name <- readRDS(sprintf("../input/Gustafsdottir/%s/cmpds.names.rds", base.dir))
  Pf.cmpd$data <- plyr::join(Pf.cmpd$data, cmpd.name, by = "Image_Metadata_BROAD_ID")
  Pf.cmpd$factor_cols <- c(Pf.cmpd$factor_cols, "Image_Metadata_SOURCE_COMPOUND_NAME")
  med <- apply(Pf.cmpd$data[,Pf.cmpd$feat_cols], 2, median)
  md <- apply(Pf.cmpd$data[,Pf.cmpd$feat_cols], 2, mad)
  Pf.cmpd$data[,Pf.cmpd$feat_cols] <- scale(Pf.cmpd$data[,Pf.cmpd$feat_cols], med, md)
}

cmpd.list <- lapply(Pf.cmpd$data$Image_Metadata_BROAD_ID, function(x) sprintf("%s-%s", str_split(x, "-")[[1]][1], str_split(x, "-")[[1]][2])) %>% unlist
cr.cmpd <- cor(Pf.cmpd$data[,Pf.cmpd$feat_cols] %>% t)
rownames(cr.cmpd) <- cmpd.list
colnames(cr.cmpd) <- cmpd.list

Pf.repurp <- readRDS("../input/repurp/2016_04_01_a549_48hr_batch1_normalized.rds")
moa.list <- Pf.repurp %>% dplyr::filter(Metadata_pert_id %in% cmpd.list) %>% dplyr::select(one_of(c("Metadata_pert_id", "Metadata_moa"))) %>% unique
v <- moa.list$Metadata_moa
same.moa <- outer(v, v, "==")
rownames(same.moa) <- moa.list$Metadata_pert_id
colnames(same.moa) <- moa.list$Metadata_pert_id
same.moa.list <- same.moa %>% melt %>% dplyr::filter(as.character(Var1) < as.character(Var2))
cr.cmpd.res <- cr.cmpd[moa.list$Metadata_pert_id, moa.list$Metadata_pert_id]
cr.cmpd.res <- cr.cmpd.res %>% melt %>% dplyr::filter(as.character(Var1) < as.character(Var2))
set0 <- cr.cmpd.res$value[which(!same.moa.list$value)]
set1 <- cr.cmpd.res$value[which(same.moa.list$value)]
t.test(set0, set1, "less")

thr <- quantile(cr.cmpd.res$value, 0.95)
v11 <- which((cr.cmpd.res$value) > thr & same.moa.list$value) %>% length
v12 <- which((cr.cmpd.res$value) > thr & !same.moa.list$value) %>% length
v21 <- which((cr.cmpd.res$value) < thr & same.moa.list$value) %>% length
v22 <- which((cr.cmpd.res$value) < thr & !same.moa.list$value) %>% length
V <- rbind(c(v11, v12), c(v21, v22))
print(V)
fisher.test(V, alternative = "greater")

u <- cr.cmpd.res[which((cr.cmpd.res$value) > thr & same.moa.list$value), ]
colnames(u)[1] <- "Metadata_pert_id"
u <- plyr::join(u, moa.list)

if (!cdp2) {
  brd.to.cid <- Pf.cmpd$data[,c("Image_Metadata_BROAD_ID", "Image_Metadata_SOURCE_COMPOUND_NAME")]
  colnames(brd.to.cid) <- c("Broad.ID", "Compound.Name")
} else {
  brd.to.cid <- read.csv("../input/CDP2/brd.to.cid.csv")  
}

brd.to.cid <- brd.to.cid[,c("Broad.ID", "Compound.Name")]
brd.to.cid$Broad.ID <- lapply(brd.to.cid$Broad.ID, function(x) sprintf("%s-%s", str_split(x, "-")[[1]][1], str_split(x, "-")[[1]][2])) %>% unlist

colnames(u)[1] <- "Broad.ID"
u <- plyr::join(u, brd.to.cid, by = "Broad.ID")
colnames(u)[1] <- "Compound 1"
colnames(u)[2] <- "Broad.ID"
u <- plyr::join(u, brd.to.cid, by = "Broad.ID")
colnames(u)[3] <- "Corr."
colnames(u)[2] <- "Compound 2"
colnames(u)[5] <- "Name 1"
colnames(u)[6] <- "Name 2"
u %>% dplyr::arrange(Metadata_moa) %>% dplyr::mutate(Corr. = round(Corr., 2)) %>% htmlTable::htmlTable()

##############

ag.ant <- function(x, y) {
  if (is.na(x) || is.na(y)) {
    return(F)
  }
  x.s <- str_split(x, " ")[[1]]
  y.s <- str_split(y, " ")[[1]]
  res <- F
  if (x.s[length(x.s)] == "agonist" && y.s[length(y.s)] == "antagonist" || 
      x.s[length(x.s)] == "antagonist" && y.s[length(y.s)] == "agonist") {
    if (paste(x.s[1:(length(x.s)-1)], collapse = " ") == paste(y.s[1:(length(y.s)-1)], collapse = " ")) {
      res <- T
    }
  }
  return(res)
}

v <- rep(1, length(moa.list$Metadata_pert_id))
same.moa <- outer(v, v, "*")
rownames(same.moa) <- 1:length(moa.list$Metadata_pert_id)
colnames(same.moa) <- 1:length(moa.list$Metadata_pert_id)
same.moa.list <- same.moa %>% melt %>% dplyr::filter(moa.list$Metadata_pert_id[Var1] < moa.list$Metadata_pert_id[Var2]) 
same.moa.list[,3] <- apply(same.moa.list, 1, function(x) ag.ant(moa.list$Metadata_moa[x[1]], moa.list$Metadata_moa[x[2]]))
same.moa <- outer(v, v, "*")
rownames(same.moa) <- moa.list$Metadata_pert_id
colnames(same.moa) <- moa.list$Metadata_pert_id
same.moa.list[,c(1, 2)] <- (melt(same.moa) %>% dplyr::filter(as.character(Var1) < as.character(Var2)))[,c(1, 2)]

thr <- quantile(cr.cmpd.res$value, 0.05)
v11 <- which((cr.cmpd.res$value) < thr & same.moa.list$value) %>% length
v12 <- which((cr.cmpd.res$value) < thr & !same.moa.list$value) %>% length
v21 <- which((cr.cmpd.res$value) > thr & same.moa.list$value) %>% length
v22 <- which((cr.cmpd.res$value) > thr & !same.moa.list$value) %>% length
V <- rbind(c(v11, v12), c(v21, v22))
print(V)
fisher.test(V, alternative = "greater")

u <- cr.cmpd.res[which((cr.cmpd.res$value) < thr & same.moa.list$value), ]
colnames(u)[1] <- "Metadata_pert_id"
u <- plyr::join(u, moa.list)
colnames(u)[1] <- "Metadata_pert_id 1"
colnames(u)[2] <- "Metadata_pert_id"
colnames(u)[4] <- "MOA 1"

u <- plyr::join(u, moa.list, by = "Metadata_pert_id")
colnames(u)[5] <- "MOA 2"

if (!cdp2) {
  brd.to.cid <- Pf.cmpd$data[,c("Image_Metadata_BROAD_ID", "Image_Metadata_SOURCE_COMPOUND_NAME")]
  colnames(brd.to.cid) <- c("Broad.ID", "Compound.Name")
} else {
  brd.to.cid <- read.csv("../input/CDP2/brd.to.cid.csv")  
}

brd.to.cid <- brd.to.cid[,c("Broad.ID", "Compound.Name")]
brd.to.cid$Broad.ID <- lapply(brd.to.cid$Broad.ID, function(x) sprintf("%s-%s", str_split(x, "-")[[1]][1], str_split(x, "-")[[1]][2])) %>% unlist

colnames(u)[1] <- "Broad.ID"
u <- plyr::join(u, brd.to.cid, by = "Broad.ID")
colnames(u)[1] <- "Compound 1"
colnames(u)[2] <- "Broad.ID"
u <- plyr::join(u, brd.to.cid, by = "Broad.ID")
colnames(u)[3] <- "Corr."
colnames(u)[2] <- "Compound 2"
colnames(u)[6] <- "Name 1"
colnames(u)[7] <- "Name 2"
u %>% dplyr::arrange(`MOA 1`) %>% dplyr::group_by(`Compound 1`, `Compound 2`) %>% dplyr::slice(1) %>% dplyr::ungroup() %>% dplyr::mutate(Corr. = round(Corr., 2)) %>% htmlTable::htmlTable()