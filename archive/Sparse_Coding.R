library(Rglpk)

find.targets <- function(cmpd, Pf.cmpd, Pf.gene, gamma = 0.01, thr = 0.95) {
  z <- Pf.cmpd$data[which(Pf.cmpd$data$Image_Metadata_BROAD_ID == cmpd), Pf.cmpd$feat_cols] %>% as.matrix() %>% as.vector()
  Y <- t(Pf.gene$data[, Pf.gene$feat_cols])
  
  w <- rep(1, NCOL(Y))
  n.var <- NROW(Y)
  
  Y <- cbind(Y, diag(rep(-1, NROW(Y))), diag(rep(0, NROW(Y))), matrix(0, n.var, length(w)))
  
  Y <- rbind(Y, cbind(matrix(0, n.var, length(w)), diag(rep(1, n.var)), diag(rep(-1, n.var)), matrix(0, n.var, length(w))))
  Y <- rbind(Y, cbind(matrix(0, n.var, length(w)), diag(rep(-1, n.var)), diag(rep(-1, n.var)), matrix(0, n.var, length(w))))
  
  Y <- rbind(Y, cbind(diag(rep(1, length(w))), matrix(0, length(w), n.var), matrix(0, length(w), n.var), diag(rep(-1, length(w)))))
  Y <- rbind(Y, cbind(diag(rep(-1, length(w))), matrix(0, length(w), n.var), matrix(0, length(w), n.var), diag(rep(-1, length(w)))))
  
  we <- c(rep(0, length(w)), rep(0, n.var), rep(1/n.var, n.var), rep(gamma, length(w)))
  ze <- c(z, rep(0, n.var), rep(0, n.var), rep(0, 2*length(w)))
  dire <- c(rep("==", length(z)), rep("<=", 2*n.var), rep("<=", 2*length(w)))
  
  bounds <- list(lower = list(ind = 1:length(we), val = c(rep(-Inf, length(w)), rep(-Inf, n.var), rep(0, n.var), rep(0, length(w)))), upper = list(ind = 1:length(we), val = c(rep(Inf, length(w)), rep(Inf, n.var), rep(Inf, n.var), rep(Inf, length(w)))))
  
  sol <- Rglpk_solve_LP(we, Y, dire, ze, bounds = bounds, max = F, control = list(verbose = T)) 
  w.sol <- sol$solution[1:length(w)]
  w.ord <- order(abs(w.sol), decreasing = T)
  ws <- w.sol[w.ord]
  v <- lapply(1:length(ws), function(i) (sum(ws[1:i])/sum(ws))) %>% unlist
  print(ws)
  
  if (v[1] < thr) {
    ind <- max(which(v < thr))
  } else {
    ind <- 1  
  }
  
  return(Pf.gene$data$Treatment[w.ord[1:ind]])
}

gn.pred <- lapply(Pf.cmpd$data$Image_Metadata_BROAD_ID, function(x) find.targets(x, Pf.cmpd, Pf.gene, 0.02, 0.95))
pred <- cr * 0
i <- 1
for (cmpd in Pf.cmpd$data$Image_Metadata_BROAD_ID) {
  pred[cmpd, gn.pred[[i]]] <- 1  
  i <- i + 1
}

d1 <- gene.chem.conn %>% melt %>% dplyr::filter((Var1 %>% as.character()) < (Var2 %>% as.character()))
d2 <- pred %>% melt %>% dplyr::filter((Var1 %>% as.character()) < (Var2 %>% as.character()))
colnames(d1)[3] <- "verified"
colnames(d2)[3] <- "predicted"

d <- plyr::join(d1, d2, by = c("Var1", "Var2"))

v11 <- d %>% dplyr::filter(predicted == 1 & verified == 1) %>% NROW
v12 <- d %>% dplyr::filter(predicted == 1 & verified == 0) %>% NROW
v21 <- d %>% dplyr::filter(predicted == 0 & verified == 1) %>% NROW
v22 <- d %>% dplyr::filter(predicted == 0 & verified == 0) %>% NROW
V <- rbind(c(v11, v12), c(v21, v22))
rownames(V) <- c("predicted", "remainder")
colnames(V) <- c("verified", "non-verified")
fisher.test(V, alternative = "greater")
print(V)