library(EBImage)

rescale.image <- function(fl.path) {
  x <- EBImage::readImage(fl.path)
  
  thr.1 <- quantile(x[,,1], 0.9999)
  thr.2 <- quantile(x[,,2], 0.9999)
  thr.3 <- quantile(x[,,3], 0.9999)

  x1 <- x[,,1]
  x2 <- x[,,2]
  x3 <- x[,,3]
  
  x1[x1 > thr.1] <- thr.1
  x2[x2 > thr.2] <- thr.2
  x3[x3 > thr.3] <- thr.3

  x[,,1] <- x1/thr.1
  x[,,2] <- x2/thr.2
  x[,,3] <- x3/thr.3
  
  x <- x - min(x)
  
  return(x)
}

dir.path <- "../results/manual/misc/compound"

lst <- list.dirs(dir.path)

for (lsti in lst[2:length(lst)]) {
  img.list <- list.files(lsti)
  img.list <- img.list[which(str_detect(img.list, ".tif"))]
  
  img.coll <- foreach (img.fl = img.list) %do% {
    img <- EBImage::readImage(paste0(lsti, "/", img.fl))
    img <- EBImage::channel(img, "rgb")
    q01 <- quantile(img[,,1], 0.01)
    q99 <- quantile(img[,,1], 0.99)
    img[img < q01] <- q01
    img[img > q99] <- q99
    
    img <- (img - q01)/q99
    img
  }
  
  x <- 0 * img.coll[[1]]
  
  x[,,1] <- img.coll[[2]][,,1] * 0.7 + img.coll[[3]][,,1] * 0.2 + img.coll[[4]][,,1] * 0.1 + img.coll[[5]][,,1] * 0.1
  x[,,2] <- img.coll[[2]][,,1] * 0.1 + img.coll[[3]][,,1] * 0.1 + img.coll[[4]][,,1] * 0.2 + img.coll[[5]][,,1] * 0.7
  x[,,3] <- img.coll[[1]][,,1] * 1
    
  fl.path <- paste0(lsti, "/",  "combined.png")
  
  EBImage::writeImage(x, fl.path)
}
