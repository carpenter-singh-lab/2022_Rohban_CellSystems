library(dplyr)
library(reshape2)
library(doParallel)

same.moa <- function(moa.list.1, moa.list.2) {
  if (is.na(moa.list.1) || is.na(moa.list.2) || moa.list.1 == "" || moa.list.2 == "")
    return(FALSE)
  moas.1 <- strsplit(moa.list.1, ", ")[[1]]
  moas.2 <- strsplit(moa.list.2, ", ")[[1]]
  return(any(moas.1 %in% moas.2) | any(moas.2 %in% moas.1))
}

same.moa <- Vectorize(same.moa)

enrichment_top_conn <- function(sm, metadata, top.perc = 0.95) {
  sm <- sm %>%
    reshape2::melt() %>%
    filter(as.character(Var1) < as.character(Var2) &
             Var1 != "DMSO" &
             Var2 != "DMSO") %>%
    left_join(.,
              metadata,
              by = c("Var1" = "Metadata_broad_sample")) %>%
    left_join(.,
              metadata,
              by = c("Var2" = "Metadata_broad_sample")) %>%
    mutate(same.moa = same.moa(Metadata_moa.x, Metadata_moa.y))

  thr <- quantile(sm$value, top.perc, na.rm = T)

  v11 <- sm %>%
    filter(value > thr & same.moa) %>%
    NROW

  v12 <- sm %>%
    filter(value > thr & !same.moa) %>%
    NROW

  v21 <- sm %>%
    filter(value < thr & same.moa) %>%
    NROW

  v22 <- sm %>%
    filter(value < thr & !same.moa) %>%
    NROW

  return(fisher.test(x = rbind(c(v11, v12), c(v21, v22)),
                     alternative = "greater"))
}

moa_recall <- function(sm, metadata, n.cores = 1, N = 1000) {
  doParallel::registerDoParallel(cores = n.cores)

  moas <- unlist(lapply(metadata$Metadata_moa, function(x) str_split(x, "\\|")[[1]]))
  moas <- unique(moas)
  moas <- setdiff(moas, NA)
  moas <- setdiff(moas, "")

  group_recall <- function(sm, brds) {
    #sm[brds, brds] %>%
    #  as.dist() %>%
    #  mean(., na.rm = T)
    x <- sm[brds, brds]
    diag(x) <- NA
    median(apply(x, 1, function(y) median(y, na.rm = T)), na.rm = T)
  }

  contains.moa <- function(moa.list, moa) {
    moa %in% str_split(moa.list, ", ")[[1]]
  }

  contains.moa <- Vectorize(contains.moa)
  brds.ref <- colnames(sm)

  res <- foreach (moa = moas, .combine = rbind) %dopar% {
    brds.moa <- metadata %>%
      filter(contains.moa(Metadata_moa, moa)) %>%
      select(Metadata_broad_sample) %>%
      as.matrix() %>%
      as.vector()

    brds.moa <- intersect(brds.moa, brds.ref)

    if (length(brds.moa) < 2) {
      return(data.frame(MOA = moa, p.value = NA))
    }

    x <- group_recall(sm = sm, brds = brds.moa)

    nulls <- foreach (j = 1:N, .combine = rbind) %do% {
      brds <- metadata %>%
        sample_n(NROW(brds.moa)) %>%
        select(Metadata_broad_sample) %>%
        as.matrix() %>%
        as.vector()

      brds <- intersect(brds, brds.ref)

      group_recall(sm = sm, brds = brds)
    }

    data.frame(MOA = moa, p.value = 1 - ecdf(nulls)(x))
  }

  return(res)
}

cmpd_classification <- function(sm, metadata, k0 = 5) {
  sm <- sm %>%
    reshape2::melt() %>%
    filter(Var1 != Var2 &
             Var1 != "DMSO" &
             Var2 != "DMSO") %>%
    left_join(.,
              metadata,
              by = c("Var1" = "Metadata_broad_sample")) %>%
    left_join(.,
              metadata,
              by = c("Var2" = "Metadata_broad_sample")) %>%
    mutate(same.moa = same.moa(Metadata_moa.x, Metadata_moa.y))

  cmpd.true.pos <- sm %>%
    arrange(-value) %>%
    group_by(Var1) %>%
    slice(1:k0) %>%
    summarise(pass = any(same.moa)) %>%
    ungroup() %>%
    select(pass) %>%
    as.matrix() %>%
    as.vector() %>%
    sum

  return(cmpd.true.pos/(sm$Var1 %>% unique %>% length))
}