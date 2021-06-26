feat.var <- feats
grp.var <- "Metadata_pert_id"
x <- x %>% dplyr::filter(Metadata_broad_sample_type == "trt")

fn <- function(v) {
  cr <- cor(t(v[,feat.var]))
  not.eq.well <- outer(v$Metadata_Well, v$Metadata_Well, "==")
  vl <- cr[not.eq.well]
  res <- median(vl[which(vl != 1)], na.rm = T)
  return(res)
}

fn2 <- function(v) {
  cr <- cor(t(v[,feat.var]))
  not.eq.well <- outer(v$Metadata_Well, v$Metadata_Well, "!=")
  res <- median(cr[not.eq.well], na.rm = T)
  return(res)
}

u.diff.well <- x %>% 
  dplyr::group_by_(grp.var) %>%
  do(data.frame(cr = fn2(.)))

u.same.well <- x %>% 
  dplyr::group_by_(grp.var) %>%
  do(data.frame(cr = fn(.)))


thr <- non.rep.cor(list(data = x, feat_col = feats, factor_col = metadata), "Metadata_Treatment", feats)

((u.same.well$cr > thr) %>% which %>% length)/(NROW(u.same.well))
((u.diff.well$cr > thr) %>% which %>% length)/(NROW(u.diff.well))