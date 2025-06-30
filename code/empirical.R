# library(corHMM)
setwd("~/corHMM/")
devtools::load_all()
setwd("~/corhmm-dredge/")

data(primates)
phy <- multi2di(primates[[1]])
phy$edge.length <- phy$edge.length + 1e-7
data <- primates[[2]]

# dredge_fits <- corHMMDredge(phy = phy, data = data, max.rate.cat = 1, 
#   pen.type = "l1", root.p = "maddfitz", lambda = 1, nstarts = 10, n.cores = 10)

# dredge_fitsSA <- corHMMDredgeSA(phy = phy, data = data, max.rate.cat = 2,
#   pen.type = "l1", root.p = "maddfitz", lambda = 1, nstarts = 5, n.cores = 5,
#   verbose = TRUE, max.iterations = 500, return.all = TRUE, initial.temp = 4)
# saveRDS(dredge_fitsSA, file = "empirical_results/empirical.RDS")
dredge_fitsSA <- readRDS("empirical_results/empirical.RDS")

new_init <- dredge_fitsSA$all_models[[14]]$index.mat
new_init[7,3] <- NA
# dredge_fitsSA2 <- corHMMDredgeSA(phy = phy, data = data, init.rate.cat = 2, max.rate.cat = 2,
#   pen.type = "l1", root.p = "maddfitz", lambda = 1, nstarts = 5, n.cores = 5,
#   verbose = TRUE, max.iterations = 500, return.all = TRUE, index_mat = new_init)
# saveRDS(dredge_fitsSA2, file = "empirical_results/empirical2.RDS")
dredge_fitsSA2 <- readRDS("empirical_results/empirical2.RDS")

getModelTable(dredge_fitsSA2$all_models)
dredge_fitsSA2$all_models[[8]]

# model_table <- getModelTable(dredge_fits)
getModelTable(dredge_fitsSA$all_models)
dredge_fitsSA$all_models[[14]]
pdf("plots/sa_info.pdf", height = 4, width = 12)
source("code/plot_SA_info.R")
dev.off()


model_table <- getModelTable(dredge_fitsSA2$all_models)
dredge_model <- dredge_fitsSA2$all_models[[8]]
# k_fold_res <- kFoldCrossValidation(dredge_model, k = 5, lambdas = c(0,0.25,0.5,0.75,1))
# saveRDS(k_fold_res, file = "empirical_results/empirical_kfold.RDS")
dredge_model$states[1,]
k_fold_res <- readRDS("empirical_results/empirical_kfold.RDS")
cv_table <- getCVTable(k_fold_res)

# profile_results_dredge <- get_batch_profile_lik(dredge_model, dredge = TRUE,
#   range_factor = 100, n_points = 200, verbose = TRUE, ncores = 8)
# saveRDS(profile_results_dredge, file = "empirical_results/empirical_profile.RDS")
profile_results_dredge <- readRDS("empirical_results/empirical_profile.RDS")
plot_batch_profile_lik(profile_results_dredge)
profile_results_dredge[[1]]$profile_table
profile_results_dredge[[2]]$profile_table

get_profile_confidence_intervals(profile_results_dredge)

dev.off()
RColorBrewer::display.brewer.all()
piecolors <- RColorBrewer::brewer.pal(8, "Paired")[c(2,8,6,4,1,7,5,3)]

anc_recon <- ancRECON(phy = dredge_model$phy, root.p = "maddfitz",
  data = dredge_model$data, 
  p = MatrixToPars(dredge_model), 
  method = "marginal", 
  rate.cat = 2, 
  rate.mat = dredge_model$index.mat, 
  get.tip.states = TRUE, collapse = FALSE)

phy <- ladderize(dredge_model$phy)

pdf("plots/asr_phy.pdf", height = 4, width = 6)
plot(phy, show.tip.label = FALSE, direction = "downwards")
tiplabels(pie = anc_recon$lik.tip.states, piecol = piecolors, cex = 0.35)
nodelabels(pie = anc_recon$lik.anc.states, piecol = piecolors, cex = 0.35)
legend("topleft", legend = colnames(dredge_model$states),
  pch=21, pt.bg = piecolors, cex = 0.6, bty="n", title = "")
text(x = -1, y = 50, label = "Estrus display | Multimale mating",
  cex = 0.6, adj=0)
dev.off()

# plot(phy, no.margin = TRUE, cex = 0.5)
plot(phy, show.tip.label = TRUE, cex = 0.5, direction = "downwards")
