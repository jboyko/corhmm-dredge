setwd("~/corHMM//")
source("~/corhmm-dredge/code/utils.R")
devtools::load_all()
setwd("~/corhmm-dredge/")

#### #### #### #### #### dep model
results_files <- dir("~/corhmm-dredge/structure_results/dep_model/", full.names = TRUE)
results <- lapply(results_files, readRDS)
tmp_dat <- data.frame(sp = NA, expand.grid(list(c("Ovi", "Vivi"), c("Warm", "Cold"))))
index_mat <- getStateMat4Dat(tmp_dat)$rate.mat
results <- lapply(results, function(x) 
  t(data.frame(
    l0 = model_test_dep(
      x$corhmm_fits$l0[[which.min(getModelTable(x$corhmm_fits$l0)$dAIC)]], 
      index_mat),
    l1 = model_test_dep(
      x$corhmm_fits$l1[[which.min(getModelTable(x$corhmm_fits$l1)$dAIC)]], 
      index_mat),
    l2 = model_test_dep(
      x$corhmm_fits$l2[[which.min(getModelTable(x$corhmm_fits$l2)$dAIC)]], 
      index_mat),
    er = model_test_dep(
      x$corhmm_fits$er[[which.min(getModelTable(x$corhmm_fits$er)$dAIC)]], 
      index_mat),
    sa = model_test_dep(
      x$corhmm_fits$sa[[which.min(getModelTable(x$corhmm_fits$sa)$dAIC)]], 
      index_mat)
  )))

# Combine results if needed
result_df_all <- do.call(rbind, results)

# how consistent are the results with the hypothesis?
test_summ <- t(data.frame(
  l0=colSums(result_df_all[rownames(result_df_all) == "l0",])/(dim(result_df_all)[1]/5),
  l1=colSums(result_df_all[rownames(result_df_all) == "l1",])/(dim(result_df_all)[1]/5),
  l2=colSums(result_df_all[rownames(result_df_all) == "l2",])/(dim(result_df_all)[1]/5),
  er=colSums(result_df_all[rownames(result_df_all) == "er",])/(dim(result_df_all)[1]/5),
  sa=colSums(result_df_all[rownames(result_df_all) == "sa",])/(dim(result_df_all)[1]/5)
))
test_summ
write.csv(test_summ, "tables/test_summary_dep.csv")

mod_names <- rownames(result_df_all)
result_df_all <- as.data.frame(result_df_all, row.names = NA)
result_df_all$model <- mod_names

#### #### #### #### #### ord model
# Create a list of simulation numbers
sim_nums <- 1:100

results_files <- dir("~/corhmm-dredge/structure_results/ord_model/", full.names = TRUE)
results <- lapply(results_files, readRDS)
tmp_dat <- data.frame(sp = NA, 
  Self = factor(c("outx", "fac", "self"), c("outx", "fac", "self")))
index_mat <- getStateMat4Dat(tmp_dat)$rate.mat
index_mat <- dropStateMatPars(index_mat, c(5,2,4))
results <- lapply(results, function(x) 
  t(data.frame(
    l0 = model_test_ord(
      x$corhmm_fits$l0[[which.min(getModelTable(x$corhmm_fits$l0)$dAIC)]], 
      index_mat),
    l1 = model_test_ord(
      x$corhmm_fits$l1[[which.min(getModelTable(x$corhmm_fits$l1)$dAIC)]], 
      index_mat),
    l2 = model_test_ord(
      x$corhmm_fits$l2[[which.min(getModelTable(x$corhmm_fits$l2)$dAIC)]], 
      index_mat),
    er = model_test_ord(
      x$corhmm_fits$er[[which.min(getModelTable(x$corhmm_fits$er)$dAIC)]], 
      index_mat),
    sa = model_test_ord(
      x$corhmm_fits$sa[[which.min(getModelTable(x$corhmm_fits$sa)$dAIC)]], 
      index_mat)
  )))

# Combine results if needed
result_df_all <- do.call(rbind, results)

# how consistent are the results with the hypothesis?
test_summ <- data.frame(
  l0=colSums(result_df_all[rownames(result_df_all) == "l0",])/(dim(result_df_all)[1]/5),
  l1=colSums(result_df_all[rownames(result_df_all) == "l1",])/(dim(result_df_all)[1]/5),
  l2=colSums(result_df_all[rownames(result_df_all) == "l2",])/(dim(result_df_all)[1]/5),
  er=colSums(result_df_all[rownames(result_df_all) == "er",])/(dim(result_df_all)[1]/5),
  sa=colSums(result_df_all[rownames(result_df_all) == "sa",])/(dim(result_df_all)[1]/5)
)
test_summ
write.csv(test_summ, "tables/test_summary_ord.csv")


#### #### #### #### #### hmm model
sim_nums <- 1:100

results_files <- dir("~/corhmm-dredge/structure_results/hmm_model/", full.names = TRUE)
results <- lapply(results_files, readRDS)
tmp_dat <- data.frame(sp = NA, expand.grid(list(c("muted", "vivid"), c("0", "1"))))
index_mat <- getStateMat4Dat(tmp_dat)$rate.mat
index_mat <- equateStateMatPars(index_mat, list(c(2,4,5,7)))
results <- lapply(results, function(x) 
  t(data.frame(
    l0 = model_test_hmm(
      x$corhmm_fits$l0[[which.min(getModelTable(x$corhmm_fits$l0)$dAIC)]], 
      index_mat),
    l1 = model_test_hmm(
      x$corhmm_fits$l1[[which.min(getModelTable(x$corhmm_fits$l1)$dAIC)]], 
      index_mat),
    l2 = model_test_hmm(
      x$corhmm_fits$l2[[which.min(getModelTable(x$corhmm_fits$l2)$dAIC)]], 
      index_mat),
    er = model_test_hmm(
      x$corhmm_fits$er[[which.min(getModelTable(x$corhmm_fits$er)$dAIC)]], 
      index_mat),
    sa = model_test_hmm(
      x$corhmm_fits$sa[[which.min(getModelTable(x$corhmm_fits$sa)$dAIC)]], 
      index_mat)
  )))

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
test_summ
write.csv(test_summ, "tables/test_summary_hmm.csv")
round(test_summ, 3) * 100
