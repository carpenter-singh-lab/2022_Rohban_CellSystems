for (i in c("_RNA", "_AGP", "_DNA", "_Mito", "_ER")) {
	ft.keyword <- i
	source("moa_enrichment_test.R", local = FALSE)
}
