# library(corHMM)
library(geiger)
library(MASS)

setwd("~/corHMM//")
source("~/corhmm-dredge/code/utils.R")
devtools::load_all()

run_simulation <- function(sim_num, save_results = TRUE) {
  # Simulate a phylogenetic tree
  phy <- sim.bdtree(n = 250)
  phy$edge.length <- phy$edge.length / max(branching.times(phy))
  
  # Set up index matrix
  tmp_dat <- data.frame(sp = NA, expand.grid(list(c("muted", "vivid"), c("0", "1"))))
  index_mat <- getStateMat4Dat(tmp_dat)$rate.mat
  index_mat <- equateStateMatPars(index_mat, list(c(2,4,5,7)))
  rate_mat <- index_mat
  p <- c(10, 10, .01, .01, 1)
  for (i in seq_along(p)) {
    rate_mat[index_mat == i] <- p[i]
  }
  rate_mat[is.na(rate_mat)] <- 0
  diag(rate_mat) <- -rowSums(rate_mat)
  # Simulate data
  sim_data <- get_sim_data(phy, rate_mat, index_mat, root.p = c(1, 0, 0, 0))
  cor_data <- get_formatted_data(sim_data$TipStates, index_mat)
  # table(sim_data$TipStates)
  # table(cor_data$V2)
  cor_data <- cor_data[,c(1,2)]
  # par(mar=c(.1,.1,.1,.1))
  # plot(phy, show.tip.label = FALSE, type="fan")
  # tiplabels(pch = 16, col = c("grey", "darkred")[as.numeric(cor_data$V2)], cex=0.55)
  # Run corHMM fits 
  corhmm_fit_l0 <- corHMM:::corHMMDredge(phy = phy, data = cor_data, max.rate.cat = 2, pen.type = "l1", root.p = "maddfitz", lambda = 0)
  corhmm_fit_l1 <- corHMM:::corHMMDredge(phy = phy, data = cor_data, max.rate.cat = 2, pen.type = "l1", root.p = "maddfitz", lambda = 1)
  corhmm_fit_l2 <- corHMM:::corHMMDredge(phy = phy, data = cor_data, max.rate.cat = 2, pen.type = "l2", root.p = "maddfitz", lambda = 1)
  corhmm_fit_er <- corHMM:::corHMMDredge(phy = phy, data = cor_data, max.rate.cat = 2, pen.type = "er", root.p = "maddfitz", lambda = 1)
  corhmm_fit_SA <- try(corHMM:::corHMMDredgeSA(phy = phy, data = cor_data, max.rate.cat = 2, pen.type = "l1", root.p = "maddfitz", lambda = 1))
  if(class(corhmm_fit_SA) == "try-error"){
    return(NULL)
  }

# single_fit <- corHMMDredgeBase(phy = phy, data = cor_data, rate.cat = 2, pen.type = "l2", root.p = "maddfitz", lambda = 1, grad = TRUE)
  
  # Select the best models
  best_l0 <- corhmm_fit_l0[[which.min(getModelTable(corhmm_fit_l0)$dAIC)]]
  best_l1 <- corhmm_fit_l1[[which.min(getModelTable(corhmm_fit_l1)$dAIC)]]
  best_l2 <- corhmm_fit_l2[[which.min(getModelTable(corhmm_fit_l2)$dAIC)]]
  best_er <- corhmm_fit_er[[which.min(getModelTable(corhmm_fit_er)$dAIC)]]
  best_sa <- corhmm_fit_SA[[which.min(getModelTable(corhmm_fit_SA)$dAIC)]]
  
  # Generate the result dataframe
  result_df <- t(data.frame(
    l0 = model_test_hmm(best_l0, index_mat),
    l1 = model_test_hmm(best_l1, index_mat),
    l2 = model_test_hmm(best_l2, index_mat),
    er = model_test_hmm(best_er, index_mat),
    sa = model_test_hmm(best_sa, index_mat)
  ))
  
  # Save intermediate results if specified
  if (save_results) {
    saveRDS(list(phy = phy, 
      sim_data = sim_data, 
      corhmm_fits = list(l0 = corhmm_fit_l0, l1 = corhmm_fit_l1, l2 = corhmm_fit_l2, er = corhmm_fit_er, sa = corhmm_fit_SA), 
      result_df = result_df),
      file = paste0("~/corhmm-dredge/structure_results/hmm_model/intermediate_results_sim_", sim_num, ".rds"))
  }
  
  # Return the result dataframe
  return(result_df)
}

library(parallel)
setwd("~/corhmm-dredge/")

# Determine the number of available cores
num_cores <- detectCores()

# Create a list of simulation numbers
done_index <- as.numeric(gsub("intermediate_results_sim_(\\d+)\\.rds", "\\1", 
  dir("structure_results/hmm_model/")))
sim_nums <- 1:100
sim_nums <- sim_nums[-done_index]

# Run simulations in parallel
results <- mclapply(sim_nums, run_simulation, mc.cores = 5)

results_files <- dir("~/corhmm-dredge/structure_results/hmm_model/", full.names = TRUE)
results <- lapply(results_files, readRDS)
results <- lapply(results, "[[", "result_df")

# results <- results[!unlist(lapply(results, function(x) class(x) == "try-error"))]

# Combine results if needed
result_df_all <- do.call(rbind, results)

# how consistent are the results with the hypothesis?
test_summ <- data.frame(
  l0=colSums(result_df_all[rownames(result_df_all) == "l0",], T)/(dim(result_df_all)[1]/5),
  l1=colSums(result_df_all[rownames(result_df_all) == "l1",], T)/(dim(result_df_all)[1]/5),
  l2=colSums(result_df_all[rownames(result_df_all) == "l2",], T)/(dim(result_df_all)[1]/5),
  er=colSums(result_df_all[rownames(result_df_all) == "er",], T)/(dim(result_df_all)[1]/5),
  sa=colSums(result_df_all[rownames(result_df_all) == "sa",], T)/(dim(result_df_all)[1]/5)
)
write.csv(test_summ, "tables/test_summary_hmm.csv")
round(test_summ, 3) * 100

# # dredging
# corhmm_fit_l0 <- corHMM:::corHMMDredge(phy = phy, data = cor_data, max.rate.cat = 1, pen.type = "l1", root.p="maddfitz", lambda = 0, grad=TRUE, merge.params = TRUE)
# 
# corhmm_fit_l1 <- corHMM:::corHMMDredge(phy = phy, data = cor_data, max.rate.cat = 1, pen.type = "l1", root.p="maddfitz", lambda = 1, grad=TRUE, merge.params = TRUE)
# 
# corhmm_fit_l2 <- corHMM:::corHMMDredge(phy = phy, data = cor_data, max.rate.cat = 1, pen.type = "l2", root.p="maddfitz", lambda = 1, grad=TRUE, merge.params = TRUE)
# 
# best_l0 <- corhmm_fit_l0[[which.min(getModelTable(corhmm_fit_l0)$dAIC)]]
# best_l1 <- corhmm_fit_l1[[which.min(getModelTable(corhmm_fit_l1)$dAIC)]]
# best_l2 <- corhmm_fit_l2[[which.min(getModelTable(corhmm_fit_l2)$dAIC)]]
# 
# result_df <- t(data.frame(l0=model_test_dep(best_l0, index_mat),
#   l1=model_test_dep(best_l1, index_mat),
#   l2=model_test_dep(best_l2, index_mat)))
# 
# result_df

# what about the ASR?
# getBestModels <- function(corhmm_fits){
#   best_l0 <- corhmm_fits$l0[[which.min(getModelTable(corhmm_fits$l0)$dAIC)]]
#   best_l1 <- corhmm_fits$l1[[which.min(getModelTable(corhmm_fits$l1)$dAIC)]]
#   best_l2 <- corhmm_fits$l2[[which.min(getModelTable(corhmm_fits$l2)$dAIC)]]
#   best_er <- corhmm_fits$er[[which.min(getModelTable(corhmm_fits$er)$dAIC)]]
#   best_sa <- corhmm_fits$sa[[which.min(getModelTable(corhmm_fits$sa)$dAIC)]]
#   return(list(l0=best_l0, l1=best_l1, l2=best_l2, er=best_er, sa=best_sa))
# }
# 
# setwd("~/corhmm-dredge/")
# 
# process_simulation <- function(i) {
#   focal_sim <- readRDS(paste0("structure_results/hmm_model/intermediate_results_sim_", i, ".rds"))
#   if(inherits(focal_sim, "try-error")){
#     return(NULL)
#   }
#   best_fits <- getBestModels(focal_sim$corhmm_fits)
#   tmp_dat <- data.frame(sp = NA, expand.grid(list(c("muted", "vivid"), c("0", "1"))))
#   index_mat <- getStateMat4Dat(tmp_dat)$rate.mat
#   index_mat <- equateStateMatPars(index_mat, list(c(2,4,5,7)))
#   tree <- focal_sim$phy
#   bt <- branching.times(tree)
#   true_asr_states <- factor(colnames(index_mat)[focal_sim$sim_data$NodeStates], colnames(index_mat))
#   true_asr_probs <- matrix(0, nrow=length(focal_sim$sim_data$NodeStates), 
#     ncol=dim(index_mat)[1], dimnames = list(names(focal_sim$sim_data$NodeStates), colnames(index_mat)))
#   for(j in 1:nrow(true_asr_probs)){
#     true_asr_probs[j, focal_sim$sim_data$NodeStates[j]] <- 1
#   }
#   
#   asr_dist <- matrix(0, nrow=nrow(true_asr_probs), ncol=4, 
#     dimnames = list(rownames(true_asr_probs), c("l0", "l1", "l2", "er")))
#   
#   get_hmm_div <- function(true_asr_probs, best_fit){
#     out <- numeric(dim(true_asr_probs)[1])
#     if(dim(best_fit$solution)[1]==4){
#       R1 <- mean(best_fit$solution[c(1,2), c(1,2)], na.rm = TRUE)
#       R2 <- mean(best_fit$solution[c(3,4), c(3,4)], na.rm = TRUE)
#       if(R1 > R2){
#         for(i in 1:nrow(true_asr_probs)){
#           out[i] <- js_divergence(true_asr_probs[i,], best_fit$states[i,])
#         }
#       }else{
#         for(i in 1:nrow(true_asr_probs)){
#           out[i] <- js_divergence(true_asr_probs[i,], best_fit$states[i,c(3,4,1,2)])
#         }
#       }
#     }else{
#       for(i in 1:nrow(true_asr_probs)){
#         out[i] <- js_divergence(true_asr_probs[i,], best_fit$states[i,c(1,2,1,2)]/2)
#       }
#     }
#     return(out)
#   }
#   
# 
#   asr_dist[,1] <- get_hmm_div(true_asr_probs, best_fits$l0)
#   asr_dist[,2] <- get_hmm_div(true_asr_probs, best_fits$l1)
#   asr_dist[,3] <- get_hmm_div(true_asr_probs, best_fits$l2)
#   asr_dist[,4] <- get_hmm_div(true_asr_probs, best_fits$er)
#   
#   # Fit linear models
#   lm_l0 <- lm(asr_dist[,1] ~ bt)
#   lm_l1 <- lm(asr_dist[,2] ~ bt)
#   lm_l2 <- lm(asr_dist[,3] ~ bt)
#   lm_er <- lm(asr_dist[,4] ~ bt)
#   
#   # Create the result data frame
#   asr_result_df <- data.frame(
#     label = c("l0", "l1", "l2", "er"),
#     intercept = c(coef(lm_l0)[1], coef(lm_l1)[1], coef(lm_l2)[1], coef(lm_er)[1]),
#     slope = c(coef(lm_l0)[2], coef(lm_l1)[2], coef(lm_l2)[2], coef(lm_er)[2]),
#     JS_dist = colMeans(asr_dist)
#   )
#   
#   # Return a list with the regression models and asr_result_df
#   list(asr_dist = asr_dist, lms = list(lm_l0 = lm_l0, lm_l1 = lm_l1, lm_l2 = lm_l2, lm_er = lm_er), asr_result_df = asr_result_df, bt=bt)
# }
# 
# # Apply the function to all 100 simulations using lapply
# results_list <- lapply(c(1:87, 89:97, 99:100), process_simulation)
# saveRDS(results_list, file = "summ_results/results_hmm.RDS")
# 
# # Extract all the asr_dists, lms, and asr_result_dfs
# all_asr_dists <- lapply(results_list, function(x) x$asr_dist)
# all_lms <- lapply(results_list, function(x) x$lms)
# all_asr_results <- lapply(results_list, function(x) x$asr_result_df)
# all_bts <- lapply(results_list, function(x) x$bt)
# 
# # COlors
# cols <- setNames(RColorBrewer::brewer.pal(4, "Set1"), c("l0", "l1", "l2", "er"))
# transp_cols <- sapply(cols, adjustcolor, alpha.f = 0.075)
# 
# # Plot all points and regression lines
# plot(x = NULL, y = NULL, xlim=range(unlist(all_bts)), ylim=c(0, max(unlist(lapply(all_asr_dists, max)))),
#   bty="n", ylab="JS Divergence", xlab = "Node Age")
# 
# lapply(results_list, function(res) {
#   points(res$bt, res$asr_dist[,1], col = transp_cols[1], pch=16)
#   points(res$bt, res$asr_dist[,2], col = transp_cols[2], pch=16)
#   points(res$bt, res$asr_dist[,3], col = transp_cols[3], pch=16)
#   points(res$bt, res$asr_dist[,4], col = transp_cols[4], pch=16)
#   abline(res$lms$lm_l0, col = transp_cols[1], lwd=0.5)
#   abline(res$lms$lm_l1, col = transp_cols[2], lwd=0.5)
#   abline(res$lms$lm_l2, col = transp_cols[3], lwd=0.5)
#   abline(res$lms$lm_er, col = transp_cols[4], lwd=0.5)
# })
# 
# # Calculate and plot the average regression lines
# avg_intercept_l0 <- median(sapply(all_lms, function(x) coef(x$lm_l0)[1]))
# avg_slope_l0 <- median(sapply(all_lms, function(x) coef(x$lm_l0)[2]))
# 
# avg_intercept_l1 <- median(sapply(all_lms, function(x) coef(x$lm_l1)[1]))
# avg_slope_l1 <- median(sapply(all_lms, function(x) coef(x$lm_l1)[2]))
# 
# avg_intercept_l2 <- median(sapply(all_lms, function(x) coef(x$lm_l2)[1]))
# avg_slope_l2 <- median(sapply(all_lms, function(x) coef(x$lm_l2)[2]))
# 
# avg_intercept_er <- median(sapply(all_lms, function(x) coef(x$lm_er)[1]))
# avg_slope_er <- median(sapply(all_lms, function(x) coef(x$lm_er)[2]))
# 
# abline(a = avg_intercept_l0, b = avg_slope_l0, col = cols[1], lwd=2)
# abline(a = avg_intercept_l1, b = avg_slope_l1, col = cols[2], lwd=2)
# abline(a = avg_intercept_l2, b = avg_slope_l2, col = cols[3], lwd=2)
# abline(a = avg_intercept_er, b = avg_slope_er, col = cols[4], lwd=2)
# 
# legend("topright", legend = names(cols), pch=21, pt.bg=cols)
# 
# # Combine all asr_result_df into one data frame
# combined_results <- do.call(rbind, all_asr_results)
# 
# # Calculate the average for intercepts, slopes, and JS_dist for each label
# average_results <- aggregate(cbind(intercept, slope, JS_dist) ~ label, data=combined_results, FUN=median)
# 
# # Print the average results
# print(average_results)
# average_results <- t(as.data.frame(average_results))
# 
# write.csv(average_results, "tables/structure_test_hmm_median.csv", row.names = T)
