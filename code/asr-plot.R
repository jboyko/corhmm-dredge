# COlors
cols <- setNames(RColorBrewer::brewer.pal(4, "Set1"), c("Mk", "l1", "l2", "er"))
transp_cols <- sapply(cols, adjustcolor, alpha.f = 0.075)

pdf("plots/JS-div-asr.pdf")
mar_st <- par()$mar
par(mfrow=c(3,1), mar = c(5.1,4.1,2.1,2.1))
#################################################################
############## CORR
#################################################################
results_list <- readRDS("summ_results/results_corr.RDS")
# Extract all the asr_dists, lms, and asr_result_dfs
all_asr_dists <- lapply(results_list, function(x) x$asr_dist)
all_lms <- lapply(results_list, function(x) x$lms)
all_asr_results <- lapply(results_list, function(x) x$asr_result_df)
all_bts <- lapply(results_list, function(x) x$bt)

# Plot all points and regression lines
plot(x = NULL, y = NULL, xlim=range(unlist(all_bts)), ylim=c(0, max(unlist(lapply(all_asr_dists, max)))),
  bty="n", ylab="JS Divergence", xlab = "")
lapply(results_list, function(res) {
  points(res$bt, res$asr_dist[,1], col = transp_cols[1], pch=16)
  points(res$bt, res$asr_dist[,2], col = transp_cols[2], pch=16)
  points(res$bt, res$asr_dist[,3], col = transp_cols[3], pch=16)
  points(res$bt, res$asr_dist[,4], col = transp_cols[4], pch=16)
  abline(res$lms$lm_l0, col = transp_cols[1], lwd=0.5)
  abline(res$lms$lm_l1, col = transp_cols[2], lwd=0.5)
  abline(res$lms$lm_l2, col = transp_cols[3], lwd=0.5)
  abline(res$lms$lm_er, col = transp_cols[4], lwd=0.5)
})

# Calculate and plot the average regression lines
avg_intercept_l0 <- median(sapply(all_lms, function(x) coef(x$lm_l0)[1]))
avg_slope_l0 <- median(sapply(all_lms, function(x) coef(x$lm_l0)[2]))

avg_intercept_l1 <- median(sapply(all_lms, function(x) coef(x$lm_l1)[1]))
avg_slope_l1 <- median(sapply(all_lms, function(x) coef(x$lm_l1)[2]))

avg_intercept_l2 <- median(sapply(all_lms, function(x) coef(x$lm_l2)[1]))
avg_slope_l2 <- median(sapply(all_lms, function(x) coef(x$lm_l2)[2]))

avg_intercept_er <- median(sapply(all_lms, function(x) coef(x$lm_er)[1]))
avg_slope_er <- median(sapply(all_lms, function(x) coef(x$lm_er)[2]))

abline(a = avg_intercept_l0, b = avg_slope_l0, col = cols[1], lwd=2)
abline(a = avg_intercept_l1, b = avg_slope_l1, col = cols[2], lwd=2)
abline(a = avg_intercept_l2, b = avg_slope_l2, col = cols[3], lwd=2)
abline(a = avg_intercept_er, b = avg_slope_er, col = cols[4], lwd=2)

H <- max(unlist(lapply(all_asr_dists, max)))
legend(x = 0, y = H, legend = names(cols), pch=21, pt.bg=cols, bg = "white")
text(x = -0.1, y = H*1.15, label = "a) Dependent generating model", xpd = TRUE, adj = 0)


#################################################################
############## ORD
#################################################################
results_list <- readRDS("summ_results/results_ord.RDS")
# Extract all the asr_dists, lms, and asr_result_dfs
all_asr_dists <- lapply(results_list, function(x) x$asr_dist)
all_lms <- lapply(results_list, function(x) x$lms)
all_asr_results <- lapply(results_list, function(x) x$asr_result_df)
all_bts <- lapply(results_list, function(x) x$bt)

# Plot all points and regression lines
plot(x = NULL, y = NULL, xlim=range(unlist(all_bts)), ylim=c(0, max(unlist(lapply(all_asr_dists, max)))),
  bty="n", ylab="JS Divergence", xlab = "")
lapply(results_list, function(res) {
  points(res$bt, res$asr_dist[,1], col = transp_cols[1], pch=16)
  points(res$bt, res$asr_dist[,2], col = transp_cols[2], pch=16)
  points(res$bt, res$asr_dist[,3], col = transp_cols[3], pch=16)
  points(res$bt, res$asr_dist[,4], col = transp_cols[4], pch=16)
  abline(res$lms$lm_l0, col = transp_cols[1], lwd=0.5)
  abline(res$lms$lm_l1, col = transp_cols[2], lwd=0.5)
  abline(res$lms$lm_l2, col = transp_cols[3], lwd=0.5)
  abline(res$lms$lm_er, col = transp_cols[4], lwd=0.5)
})

# Calculate and plot the average regression lines
avg_intercept_l0 <- median(sapply(all_lms, function(x) coef(x$lm_l0)[1]))
avg_slope_l0 <- median(sapply(all_lms, function(x) coef(x$lm_l0)[2]))

avg_intercept_l1 <- median(sapply(all_lms, function(x) coef(x$lm_l1)[1]))
avg_slope_l1 <- median(sapply(all_lms, function(x) coef(x$lm_l1)[2]))

avg_intercept_l2 <- median(sapply(all_lms, function(x) coef(x$lm_l2)[1]))
avg_slope_l2 <- median(sapply(all_lms, function(x) coef(x$lm_l2)[2]))

avg_intercept_er <- median(sapply(all_lms, function(x) coef(x$lm_er)[1]))
avg_slope_er <- median(sapply(all_lms, function(x) coef(x$lm_er)[2]))

abline(a = avg_intercept_l0, b = avg_slope_l0, col = cols[1], lwd=2)
abline(a = avg_intercept_l1, b = avg_slope_l1, col = cols[2], lwd=2)
abline(a = avg_intercept_l2, b = avg_slope_l2, col = cols[3], lwd=2)
abline(a = avg_intercept_er, b = avg_slope_er, col = cols[4], lwd=2)

H <- max(unlist(lapply(all_asr_dists, max)))
# legend(x = 0, y = H, legend = names(cols), pch=21, pt.bg=cols, bg = "white")
text(x = -0.1, y = H*1.15, label = "b) Ordered generating model", xpd = TRUE, adj = 0)

#################################################################
############## hmm
#################################################################
results_list <- readRDS("summ_results/results_hmm.RDS")
# Extract all the asr_dists, lms, and asr_result_dfs
all_asr_dists <- lapply(results_list, function(x) x$asr_dist)
all_lms <- lapply(results_list, function(x) x$lms)
all_asr_results <- lapply(results_list, function(x) x$asr_result_df)
all_bts <- lapply(results_list, function(x) x$bt)

# Plot all points and regression lines
plot(x = NULL, y = NULL, xlim=range(unlist(all_bts)), ylim=c(0, max(unlist(lapply(all_asr_dists, max)))),
  bty="n", ylab="JS Divergence", xlab = "Node Age")
lapply(results_list, function(res) {
  points(res$bt, res$asr_dist[,1], col = transp_cols[1], pch=16)
  points(res$bt, res$asr_dist[,2], col = transp_cols[2], pch=16)
  points(res$bt, res$asr_dist[,3], col = transp_cols[3], pch=16)
  points(res$bt, res$asr_dist[,4], col = transp_cols[4], pch=16)
  abline(res$lms$lm_l0, col = transp_cols[1], lwd=0.5)
  abline(res$lms$lm_l1, col = transp_cols[2], lwd=0.5)
  abline(res$lms$lm_l2, col = transp_cols[3], lwd=0.5)
  abline(res$lms$lm_er, col = transp_cols[4], lwd=0.5)
})

# Calculate and plot the average regression lines
avg_intercept_l0 <- median(sapply(all_lms, function(x) coef(x$lm_l0)[1]))
avg_slope_l0 <- median(sapply(all_lms, function(x) coef(x$lm_l0)[2]))

avg_intercept_l1 <- median(sapply(all_lms, function(x) coef(x$lm_l1)[1]))
avg_slope_l1 <- median(sapply(all_lms, function(x) coef(x$lm_l1)[2]))

avg_intercept_l2 <- median(sapply(all_lms, function(x) coef(x$lm_l2)[1]))
avg_slope_l2 <- median(sapply(all_lms, function(x) coef(x$lm_l2)[2]))

avg_intercept_er <- median(sapply(all_lms, function(x) coef(x$lm_er)[1]))
avg_slope_er <- median(sapply(all_lms, function(x) coef(x$lm_er)[2]))

abline(a = avg_intercept_l0, b = avg_slope_l0, col = cols[1], lwd=2)
abline(a = avg_intercept_l1, b = avg_slope_l1, col = cols[2], lwd=2)
abline(a = avg_intercept_l2, b = avg_slope_l2, col = cols[3], lwd=2)
abline(a = avg_intercept_er, b = avg_slope_er, col = cols[4], lwd=2)

H <- max(unlist(lapply(all_asr_dists, max)))
# legend(x = 0, y = H, legend = names(cols), pch=21, pt.bg=cols, bg = "white")
text(x = -0.1, y = H*1.15, label = "c) Hidden Markov generating model", xpd = TRUE, adj = 0)


# 
dev.off()





