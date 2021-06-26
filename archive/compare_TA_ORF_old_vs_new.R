load("../input/TA/Initial_analysis_workspace_old.RData")
Px.old <- Px.org.org.ta
load("../input/TA/Initial_analysis_workspace.RData")
Px.new <- Px.org.org.ta

d.old <- Px.old[, c("Plate_", "Well", "Cells_Correlation_Correlation_DNA_ER")]
d.new <- Px.new[, c("Plate", "Well", "Cells_Correlation_Correlation_DNA_ER")]
colnames(d.old)[1] <- "Plate"

d.jn <- plyr::join(d.old, d.new, by = c("Plate", "Well"))
plot(d.jn[,3], d.jn[,4])