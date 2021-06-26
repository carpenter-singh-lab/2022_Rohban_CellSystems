D <- data.frame(thr = c(), enr.ratio = c())
sq <- seq(from = 0.980, to = 0.999, by = 0.002)
pb <- progress::progress_bar$new(total = length(sq))

for (cor.thr.per.i in sq) {
	cor.thr.per <- cor.thr.per.i
	knitr::knit2html("CDRP_per_gene_analysis.Rmd", envir=parent.frame())
	enr <- readRDS("tmp.rds")
	D <- rbind(D, data.frame(thr = cor.thr.per, enr.ratio = enr))
	pb$tick()
}
saveRDS(D, "data_eval_prof.rds")
