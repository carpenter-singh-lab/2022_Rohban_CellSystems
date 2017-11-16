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

dir.path <- "../results/manual/toxic_phenotypes"

lst <- list.dirs(dir.path)

for (lsti in lst[2:length(lst)]) {
  fl.path <- paste0(lsti, "/",  "combined.png")
  
  x <- rescale.image(fl.path)
  
  EBImage::writeImage(x, fl.path)
}
