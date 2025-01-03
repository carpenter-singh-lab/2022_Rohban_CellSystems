cmpd <- "peucedanin"

#######

x <- Pf.cdrp$data %>% dplyr::filter(CPD_NAME == cmpd) %>% ungroup() %>% dplyr::select(one_of(f)) %>% slice(1)

y <- Pf.gust$data %>% dplyr::filter(str_to_lower(Image_Metadata_SOURCE_COMPOUND_NAME) == "curcumin") %>% dplyr::select(one_of(f)) %>% slice(1)

D <- data.frame(x = x %>% as.matrix() %>% as.vector(), y = y %>% as.matrix() %>% as.vector())

rownames(D) <- f

#######

brd.cdrp <- Pf.cdrp$data %>% dplyr::filter(CPD_NAME == cmpd) %>% ungroup() %>% slice(1) %>% dplyr::select(Metadata_broad_sample) %>% as.matrix() %>% as.vector()

Pf <- readRDS("../results/master/2017-06-02_4df9ed93/Pf_bio_new.rds")
u.cdrp <- rep.cor(list(data = Pf),
                  grp.var = "Metadata_broad_sample",
                  feat.var = f)


pl.well.cdrp <- Pf %>% dplyr::filter(Metadata_broad_sample == brd.cdrp) %>% dplyr::select(Metadata_Plate, Metadata_Well)

brd.gust <- Pf.gust$data %>% dplyr::filter(str_to_lower(Image_Metadata_SOURCE_COMPOUND_NAME) == cmpd) %>% ungroup() %>% slice(1) %>% dplyr::select(Image_Metadata_BROAD_ID) %>% as.matrix() %>% as.vector()

load("../input/Gustafsdottir/Initial_analysis.RData")
pl.well.gust <- Pf.full.plate.norm$data %>% dplyr::filter(Image_Metadata_BROAD_ID == brd.gust) %>% dplyr::select(Plate, Well)

cdrp.rep <- u.cdrp %>% dplyr::filter(Metadata_broad_sample == brd.cdrp)

u.gust <- rep.cor(Pf.full.plate.norm,
             grp.var = "Image_Metadata_BROAD_ID",
             feat.var = Pf.full.plate.norm$feat_cols)

gust.rep <- u.gust %>% dplyr::filter(Image_Metadata_BROAD_ID == brd.gust)
