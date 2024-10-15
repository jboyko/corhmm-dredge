# library(corHMM)
setwd("~/corHMM/")
devtools::load_all()
setwd("~/corhmm-dredge/")

data(primates)
phy <- multi2di(primates[[1]])
phy$edge.length <- phy$edge.length + 1e-7
data <- primates[[2]]
dredge_fits <- corHMMDredge(phy = phy, data = data, max.rate.cat = 1, pen.type = "l1", root.p = "maddfitz", lambda = 1, nstarts = 10, n.cores = 10)
model_table <- getModelTable(dredge_fits)

dredge_model <- dredge_fits[[which.min(model_table$dAIC)]]
k_fold_res <- kFoldCrossValidation(dredge_model, k = 5, lambdas = c(0,0.25,0.5,0.75,1))
cv_table <- getCVTable(k_fold_res)

profile_results_dredge <- get_batch_profile_lik(dredge_model, dredge = TRUE,
  range_factor = 100, n_points = 50, verbose = TRUE, ncores = 10)
plot_batch_profile_lik(profile_results_dredge)

corhmm_model <- corHMM(phy, data, 1, model = "ARD", root.p = "maddfitz", collapse = FALSE, nstarts = 10, n.cores = 10)
corhmm_model$index.mat <- dropStateMatPars(corhmm_model$index.mat, c(1,3,4,8))
profile_results_corhmm <- get_batch_profile_lik(corhmm_model, dredge = FALSE,
  range_factor = 100, n_points = 50, verbose = TRUE, ncores = 10)
plot_batch_profile_lik(profile_results_corhmm)


corhmm_model_2 <- corHMM(phy, data, 1, rate.mat = dredge_model$index.mat, root.p = "maddfitz", collapse = FALSE, nstarts = 10, n.cores = 10)
profile_results_corhmm_2 <- get_batch_profile_lik(corhmm_model_2, dredge = FALSE,
  range_factor = 100, n_points = 50, verbose = TRUE, ncores = 10)
plot_batch_profile_lik(profile_results_corhmm_2)




corhmm_model_2$AIC - corhmm_model$AIC

piecolors=c("white","black","red","yellow","forestgreen","blue","coral","aquamarine","darkorchid","gold","grey","yellow","#3288BD","#E31A1C")

dev.off()
RColorBrewer::display.brewer.all()
piecolors <- RColorBrewer::brewer.pal(4, "Paired")

par(mar=c(.1,.1,.1,.1), mfrow=c(1,2))
plot(dredge_model$phy, show.tip.label = FALSE)
tiplabels(pie = dredge_model$tip.states, cex=0.5, piecol = piecolors)
nodelabels(pie = dredge_model$states, piecol = piecolors)
legend("topleft", legend = colnames(dredge_model$states), 
  pch=21, pt.bg = piecolors, cex = 0.75, bty="n", title = "")
text(x = -1, y = 61.5, label = "Estrus display | Multimale mating", 
  cex = 0.75, adj=0)

plot(corhmm_model$phy, show.tip.label = FALSE, direction = "leftwards")
tiplabels(pie = corhmm_model$tip.states, cex=0.5, piecol = piecolors)
nodelabels(pie = corhmm_model$states, piecol = piecolors)




p <- corhmm_model$solution[!is.na(corhmm_model$index.mat)]
p[4] <- p[4]*1000
corhmm_model_3 <- corHMM(phy, data, 1, root.p = "maddfitz", collapse = FALSE, rate.mat = corhmm_model$index.mat, p = p)

par(mar=c(.1,.1,.1,.1), mfrow=c(1,2))
plot(corhmm_model_3$phy, show.tip.label = FALSE)
tiplabels(pie = corhmm_model_3$tip.states, cex=0.5, piecol = piecolors)
nodelabels(pie = corhmm_model_3$states, piecol = piecolors)
legend("topleft", legend = colnames(corhmm_model_3$states), 
  pch=21, pt.bg = piecolors, cex = 0.75, bty="n", title = "")
text(x = -1, y = 61.5, label = "Estrus display | Multimale mating", 
  cex = 0.75, adj=0)

plot(corhmm_model$phy, show.tip.label = FALSE, direction = "leftwards")
tiplabels(pie = corhmm_model$tip.states, cex=0.5, piecol = piecolors)
nodelabels(pie = corhmm_model$states, piecol = piecolors)
