library(dplyr)
library(magrittr)
library(stringr)
library(foreach)
library(doParallel)

read.dataset.single.cell.dmso <- function(data.set.name) {
  if (data.set.name == "CDRP") {
    Px <- readRDS("../results/master/2017-07-12_f947a2d1/CDRP_dmso_single_cell.rds")
    Pf <- list(data = Px, feat_cols = colnames(Px), factor_cols = NULL)
  } else if (data.set.name == "BBBC022") {
    #Px <- readRDS("../results/master/2017-07-12_2f6c36eb/BBBC022_dmso_single_cell.rds")
    Pf.gust <- readRDS("../input/Gustafsdottir/BBBC022_dmso/neg_ctrl_samples.rds")
    f2 <- colnames(Pf.gust$data)
    f2 <- str_replace_all(f2, "Syto", "RNA")
    f2 <- str_replace_all(f2, "Hoechst", "DNA")
    f2 <- str_replace_all(f2, "Ph_golgi", "AGP")
    f2 <- str_replace_all(f2, "_3", "_3_0")
    f2 <- str_replace_all(f2, "_5", "_5_0")
    colnames(Pf.gust$data) <- f2
    Pf.gust$feat_cols <- setdiff(colnames(Pf.gust$data), Pf.gust$factor_cols)
    Pf <- Pf.gust
  } else if (data.set.name == "Repurposing") {
    Px.1 <- readRDS("../results/master/2017-07-13_e204cf3a__0/Repurposing_dmso_single_cell.rds")
    Px.2 <- readRDS("../results/master/2017-07-13_e204cf3a__1/Repurposing_dmso_single_cell.rds")
    Px.3 <- readRDS("../results/master/2017-07-17_e204cf3a/Repurposing_dmso_single_cell.rds")
    ft <- intersect(intersect(colnames(Px.1), colnames(Px.2)), colnames(Px.3))
    Px <- rbind(Px.1[,ft], Px.2[,ft], Px.3[,ft])
    Pf <- list(data = Px, feat_cols = colnames(Px), factor_cols = NULL)
  }
  
  return(Pf)
}
  
read.dataset <- function(data.set.name, just.bioactives = T, dose.closest = 10, standardize.well = T) {
  if (data.set.name == "BBBC022") {
    load("../input/Gustafsdottir/Initial_analysis.RData")
    Pf.full.plate.norm$data %<>% filter(Image_Metadata_BROAD_ID != "" | 
                                          (!str_detect(Well, "01") &
                                             !str_detect(Well, "24") &
                                             !str_detect(Well, "A") &
                                             !str_detect(Well, "P")))
    metad <- readr::read_csv("../input/Gustafsdottir/BBBC022_v1_image.csv")
    metad <- metad[,c("Image_Metadata_CPD_PLATE_MAP_NAME", "Image_Metadata_PlateID")] %>% unique
    Pf.full.plate.norm$data %<>% left_join(., metad, by = c("Plate" = "Image_Metadata_PlateID"))
    Pf.full.plate.norm$factor_cols <- c(Pf.full.plate.norm$factor_cols, 
                                        "Image_Metadata_CPD_PLATE_MAP_NAME")
    
    Pf.full.plate.norm$data %<>% dplyr::mutate(Metadata_plate_well = paste(Image_Metadata_CPD_PLATE_MAP_NAME, Well, sep = "_"))
    Pf.full.plate.norm$factor_cols <- c(Pf.full.plate.norm$factor_cols, 
                                        "Metadata_plate_well")
    
    
    MOAs <- read.csv("../input/Gustafsdottir/MOAs.csv")
    Pf.full.plate.norm$data %<>% 
      dplyr::left_join(., 
                       MOAs, 
                       by = c("Image_Metadata_SOURCE_COMPOUND_NAME" = "Name"))
    
    Pf.full.plate.norm$factor_cols <- c(Pf.full.plate.norm$factor_cols, 
                                        colnames(MOAs)) %>% unique()
    
    Pf.gust <- Pf.full.plate.norm
    Pf.gust$data %<>% 
      dplyr::rename(Metadata_broad_sample = Image_Metadata_BROAD_ID,
                    Metadata_moa = MOA,
                    Metadata_pert_iname = Image_Metadata_SOURCE_COMPOUND_NAME,
                    Metadata_target = Target,
                    Metadata_plate_map_name = Image_Metadata_CPD_PLATE_MAP_NAME,
                    Metadata_Plate = Plate,
                    Metadata_Well = Well)
    
    Pf.gust$data %<>% 
      dplyr::mutate(Metadata_pert_iname = ifelse(Metadata_pert_iname == "", 
                                                 "DMSO",
                                                 as.character(Metadata_pert_iname))) %>%
      dplyr::mutate(Metadata_pert_iname = str_to_lower(Metadata_pert_iname)) 
    
    if (standardize.well) {
  	Pf.gust$data %<>% dplyr::mutate(Metadata_Well = str_to_lower(Metadata_Well))
    }
    
    Pf.gust$factor_cols <- c("Metadata_broad_sample",
                             "Metadata_moa",
                             "Metadata_pert_iname",
                             "Metadata_target",
                             "Metadata_plate_map_name",
                             "Metadata_Plate",
                             "Metadata_Well",
                             "Metadata_plate_well")
    
    f2 <- colnames(Pf.gust$data)
    f2 <- str_replace_all(f2, "Syto", "RNA")
    f2 <- str_replace_all(f2, "Hoechst", "DNA")
    f2 <- str_replace_all(f2, "Ph_golgi", "AGP")
    f2 <- str_replace_all(f2, "_3", "_3_0")
    f2 <- str_replace_all(f2, "_5", "_5_0")
    colnames(Pf.gust$data) <- f2
    Pf.gust$feat_cols <- setdiff(colnames(Pf.gust$data), Pf.gust$factor_cols)
    
    return(Pf.gust)
      
  } else if (data.set.name == "CDRP") {
    if (just.bioactives) {
      Pf <- readRDS("../results/master/2017-06-26_5b0bbf7c/Pf_bio_new.rds")
    } else {
      Pf.1 <- readRDS("../results/master/2017-06-26_5b0bbf7c/Pf_bio_new.rds")
      Pf.2 <- readRDS("../results/master/2017-06-26_5b0bbf7c/Pf_DOS_new.rds")
      Pf <- rbind(Pf.1, Pf.2)
    }    
    
    Pf %<>% dplyr::filter(Metadata_broad_sample != "DMSO" | 
                     (!str_detect(Metadata_Well, "01") &
                        !str_detect(Metadata_Well, "24") &
                        !str_detect(Metadata_Well, "a") &
                        !str_detect(Metadata_Well, "p")))
    
    Pf %<>% dplyr::mutate(Metadata_plate_well = paste(Metadata_Plate_Map_Name, Metadata_Well, sep = "_"))
    
    feat <- Pf %>% 
      dplyr::select(-dplyr::contains("Metadata_")) %>% 
      colnames()
    
    metadata.ext <- readr::read_csv("../input/CDP2/cdrp.cpd.meta.csv")
    brd.full <- unique(Pf$Metadata_broad_sample)
    brds <- lapply(brd.full, function(x) paste(str_split(x,
                                                         "-")[[1]][1:2], 
                                               sep = "-", collapse = "-")) %>% unlist()
    brd.mapping <- data.frame(BROAD_CPD_ID = brds, 
                              Image_Metadata_BROAD_ID = brd.full)
    metadata.ext %<>% 
      dplyr::select(BROAD_CPD_ID, CPD_NAME) %>% 
      dplyr::filter(BROAD_CPD_ID %in% brds) %>%
      dplyr::mutate(CPD_NAME = ifelse(str_detect(CPD_NAME, "BRD-"), "", CPD_NAME)) %>%
      dplyr::inner_join(., brd.mapping, by = "BROAD_CPD_ID") %>%
      dplyr::select(-BROAD_CPD_ID)
    
    metadata <- data.frame(Metadata_broad_sample = unique(Pf$Metadata_broad_sample))
    metadata %<>% 
      dplyr::left_join(., metadata.ext, by = c("Metadata_broad_sample" = "Image_Metadata_BROAD_ID"))
    
    Pf %<>% dplyr::left_join(., metadata, by = "Metadata_broad_sample")
    
    Pf %<>% dplyr::mutate(CPD_NAME = str_to_lower(CPD_NAME))
    Pf %<>% dplyr::mutate(CPD_NAME_san = str_replace_all(CPD_NAME, "-", "")) %>%
      dplyr::mutate(CPD_NAME_san = str_replace_all(CPD_NAME_san, " ", "")) %>%
      dplyr::mutate(CPD_NAME_san = str_replace_all(CPD_NAME_san, "\\(", "")) %>%
      dplyr::mutate(CPD_NAME_san = str_replace_all(CPD_NAME_san, "\\)", "")) 
    
    MOA <- read.csv("../input/moas.txt", sep = "\t")
    MOA %<>% dplyr::mutate(Name = str_to_lower(Name)) %>%
      dplyr::mutate(Name_san = str_replace_all(Name, "-", "")) %>%
      dplyr::mutate(Name_san = str_replace_all(Name_san, " ", "")) %>%
      dplyr::mutate(Name_san = str_replace_all(Name_san, "\\(", "")) %>%
      dplyr::mutate(Name_san = str_replace_all(Name_san, "\\)", "")) 
    
    Pf %<>% dplyr::left_join(., MOA, by = c("CPD_NAME_san" = "Name_san"))
    Pf %<>% 
      dplyr::rename(Metadata_moa = MOA,
                    Metadata_pert_iname = CPD_NAME,
                    Metadata_target = Target,
                    Metadata_plate_map_name = Metadata_Plate_Map_Name)
    
    Pf %<>% 
      dplyr::mutate(Metadata_pert_iname = ifelse(Metadata_broad_sample == "DMSO", 
                                                 "DMSO",
                                                 as.character(Metadata_pert_iname))) %>%
      dplyr::mutate(Metadata_pert_iname = ifelse(is.na(Metadata_pert_iname) | Metadata_pert_iname == "",
                                                 Metadata_broad_sample, 
                                                 Metadata_pert_iname),
                    Metadata_pert_iname = str_to_lower(Metadata_pert_iname)
                    )
    if (standardize.well) {        
	Pf %<>% dplyr::mutate(Metadata_Well = str_to_lower(Metadata_Well))
    } 

    meta.col <- setdiff(colnames(Pf), feat)
    Pf.cdrp <- list(data = Pf, feat_cols = feat, factor_cols = meta.col)
    
    return(Pf.cdrp)
  } else if (data.set.name == "Repurposing") {
    x <- readRDS("../input/repurp/2016_04_01_a549_48hr_batch1_normalized.rds")
    
    x <- cbind(x, data.frame(Metadata_Treatment = paste(x$Metadata_pert_id, x$Metadata_mg_per_ml, sep = "@")))
    feats <- colnames(x)
    feats <- feats[which(!str_detect(feats, "Metadata"))]
    metadata <- colnames(x)
    metadata <- metadata[which(str_detect(metadata, "Metadata"))]
    
    Pf.repurp <- list(data = x, 
                      feat_cols = feats, 
                      factor_cols = metadata)
    
    Pf.repurp$data %<>% dplyr::mutate(Metadata_plate_well = paste(Metadata_plate_map_name, Metadata_Well, sep = "_"))
    
    Pf.repurp$factor_cols <- c(Pf.repurp$factor_cols, "Metadata_plate_well")
    if (!is.null(dose.closest)) {
      trts <- Pf.repurp$data %>% 
        dplyr::mutate(tmp_dose = abs(Metadata_mmoles_per_liter - dose.closest)) %>%
        dplyr::arrange(tmp_dose) %>%
        dplyr::group_by(Metadata_pert_iname) %>%
        dplyr::slice(1) %>%
        dplyr::ungroup(.) %>%
        dplyr::select(Metadata_Treatment) %>%
        as.matrix() %>%
        as.vector()
    } else {
      trts <- Pf.repurp$data %>% 
        dplyr::select(Metadata_Treatment) %>%
        as.matrix() %>%
        as.vector()
    }
    
    Pf.repurp$data %<>%
      dplyr::filter(Metadata_Treatment %in% c(trts, "NA@0")) %>%
      dplyr::mutate(Metadata_pert_iname = ifelse(is.na(Metadata_pert_iname),
                                                 "DMSO",
                                                 Metadata_pert_iname),
                    Metadata_pert_iname = str_to_lower(Metadata_pert_iname))

    if (standardize.well) {        
	Pf.repurp$data %<>% dplyr::mutate(Metadata_Well = str_to_lower(Metadata_Well))
    }
    
    Pf.repurp$data %<>% filter(Metadata_pert_iname != "dmso" | 
                     (!str_detect(Metadata_Well, "01") &
                        !str_detect(Metadata_Well, "24") &
                        !str_detect(Metadata_Well, "a") &
                        !str_detect(Metadata_Well, "p")))
    
    return(Pf.repurp)
  } else if (data.set.name == "MOA") {
    lst <- list.dirs("../input/MOA")
    lst <- lst[2:length(lst)]
    data <- foreach (l = lst, .combine = rbind) %do% {
      pl.name <- str_split(l, "/")[[1]][4]
      read.csv(sprintf("%s/%s_normalized.csv", l, pl.name))
    }
    
    w <- apply(data, 2, function(x) any(is.na(x))) %>% which
    w.nm <- names(w)
    w.nm <- w.nm[which(!str_detect(w.nm, "Metadata"))]
    w.n <- setdiff(colnames(data), w.nm)
    data <- data[, w.n]
    
    data %<>% 
      dplyr::mutate(Metadata_pert_iname = str_to_lower(Metadata_broad_sample)) %>%
      dplyr::mutate(Metadata_broad_sample = paste(Metadata_broad_sample, 
                                                  Metadata_cell_id, 
                                                  Metadata_mmoles_per_liter, 
                                                  sep = "@")) %>%
      dplyr::mutate(Metadata_plate_well = paste(Metadata_Plate_Map_Name, Metadata_Well, sep = "_"))

    if (standardize.well) {        
	data %<>% dplyr::mutate(Metadata_Well = str_to_lower(Metadata_Well))
    }
    
    y <- colnames(data)
    y <- y[which(str_detect(y, "Metadata"))]
    x <- setdiff(colnames(data), y)
    
    Pf.moa <- list(data = data, feat_cols = x, factor_cols = y)
    Pf.moa$data %<>% filter(Metadata_pert_iname != "dmso" | 
                                 (!str_detect(Metadata_Well, "01") &
                                    !str_detect(Metadata_Well, "24") &
                                    !str_detect(Metadata_Well, "a") &
                                    !str_detect(Metadata_Well, "p")))
    return(Pf.moa)
  } else {
    return(NULL)
  }
}

read.dmso.single.cell <- function(data.set.name, just.bioactives = T, dose.closest = 10, well.samples = 4, n.cores = 3) {

  doParallel::registerDoParallel(cores = n.cores)
  Pf <- read.dataset(data.set.name, just.bioactives, dose.closest, F)
  plate.well <- Pf$data %>% 
    dplyr::filter(Metadata_pert_iname == "dmso") %>%
    sample_n(well.samples) %>%
    dplyr::select(Metadata_Plate, Metadata_Well) %>%
    dplyr::arrange(Metadata_Plate)
  pl <- plate.well$Metadata_Plate %>% unique
  if (data.set.name == "CDRP") {
    base.path <- "~/efs/2015_Bray_GigaScience/workspace/backend/CDRP"
  } else if (data.set.name == "Repurposing") {
    base.path <- "~/bucket2/projects/2015_10_05_DrugRepurposing_AravindSubramanian_GolubLab_Broad/workspace/backend/2016_04_01_a549_48hr_batch1"
  } else if (data.set.name == "BBBC022") {
    base.path <- "~/efs/2016_12_13_Cytominer_Janssen/workspace/backend/BBBC022_2013"
  } else {
    return(NULL)
  }
 
  print(plate.well)
 
  prf <- foreach (pli = pl, .combine = rbind) %dopar% {
    if (data.set.name == "Repurposing") {
	system("mkdir -p tmp")
    	system(sprintf("aws s3 sync s3://imaging-platform/projects/2015_10_05_DrugRepurposing_AravindSubramanian_GolubLab_Broad/workspace/backend/2016_04_01_a549_48hr_batch1/%s tmp --exclude \"*.csv\" --exclude \"*.gct\" ", pli))
	db <- dplyr::src_sqlite(sprintf("tmp/%s.sqlite", pli))
    } else {
    	db <- dplyr::src_sqlite(sprintf("%s/%s/%s.sqlite", base.path, pli, pli))
    }
    
    cells <- dplyr::tbl(db, "Cells") 
    cytoplasm <- dplyr::tbl(db, "Cytoplasm") 
    nuclei <- dplyr::tbl(db, "Nuclei") 
    img <- dplyr::tbl(db, "Image")

    wells <- plate.well %>%
      dplyr::filter(Metadata_Plate == pli) %>%
      dplyr::select(Metadata_Well) %>%
      as.matrix() %>%
      as.vector()

    print(wells)
    prf <- foreach (well = wells, .combine = rbind) %do% {
	    prf <- img %>%
		dplyr::filter(Image_Metadata_Well == well) %>%
		dplyr::select(ImageNumber) %>%
		dplyr::left_join(., cells, by = "ImageNumber") %>%
		dplyr::left_join(., cytoplasm, by = c("ImageNumber", "ObjectNumber")) %>%
		dplyr::left_join(., nuclei, by = c("ImageNumber", "ObjectNumber")) %>%
	dplyr::collect()
prf
    }

#    wells <- plate.well %>% 
#      dplyr::filter(Metadata_Plate == pli) %>%
#      dplyr::select(Metadata_Well)


 #    img <- dplyr::tbl(db, "Image")

 #    image <- wells %>%
 #     dplyr::left_join(., img, by = c("Metadata_Well" = "Image_Metadata_Well"), copy = T) %>%
 #     dplyr::select(ImageNumber) %>%
 #     dplyr::collect() %>%
 #     as.matrix() %>%
 #     as.vector()

 #    print(image)
 #    prf <- foreach (img.number = image, .combine = rbind) %do% {
#	c1 <- cells %>% dplyr::filter(ImageNumber == img.number) 
#        c2 <- cytoplasm %>% dplyr::filter(ImageNumber == img.number) 
#	c3 <- nuclei %>% dplyr::filter(ImageNumber == img.number)

#	prf <- c1 %>%
#          dplyr::left_join(., c2, by = c("ImageNumber", "ObjectNumber")) %>%
#          dplyr::left_join(., c3, by = c("ImageNumber", "ObjectNumber")) %>%
#	  dplyr::collect()
#
#	prf
#     }

    prf
  }

  
  ft <- colnames(prf)
  ft <- ft[which(!str_detect(ft, "Metadata") & !str_detect(ft, "ImageNumber") 
       & !str_detect(ft, "ObjectNumber") & !str_detect(ft, "TableNumber"))]
  if (data.set.name == "Repurposing") {
	system("rm -R tmp")
  }
  return(prf[,ft])
}
