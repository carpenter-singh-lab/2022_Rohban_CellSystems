rm(list = ls())
source("read_dataset.R")

# processing 
dataset.name <- "Repurposing"
no.of.wells <- 3
just.bioactives <- T
dose.closest <- 10
no.cores <- 4

Px <- read.dmso.single.cell(dataset.name, just.bioactives, dose.closest, no.of.wells, no.cores)

saveRDS(Px, sprintf("%s_dmso_single_cell.rds", dataset.name))
