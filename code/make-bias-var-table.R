setwd("~/corhmm-dredge/")
library(corHMM)

source("code/utils.R")

# sim 1 
trees <- lapply(dir("trees/", full.names = TRUE), read.tree)
simulation <- "01"
index_mat <- get_index_mat(nChar=1, nStates=2, nRateClass=1)
# load
par_table_name <- paste0("par_table-", simulation, ".csv")
full_dat_name <- paste0("full_data-", simulation, ".RDS")
cor_dat_name <- paste0("cor_data-", simulation, ".RDS")
sim_res_files <- dir("res/", full.names = TRUE)[grep(paste0("res", simulation, "_"), dir("res/"))]
tmp <- get_par_table(index_mat, 100, mean = 0, sd = 0.25)  
par_table <- read.csv(paste0("parameter_tables/", par_table_name))
colnames(par_table) <- colnames(tmp)
rate_mats <- get_rate_mats(index_mat, par_table)
full_dat <- readRDS(paste0("data/", full_dat_name))
cor_dat <- lapply(full_dat, "[[", "cor_dat")
res_list <- lapply(sim_res_files, readRDS)
names(res_list) <- gsub(".*_", "", sim_res_files) %>% gsub(".RDS", "", .)
# format data and compare results
df_list <- lapply(res_list, function(x) do.call(rbind, lapply(x, get_solution_from_res)))
df_true <- do.call(rbind, lapply(full_dat, function(x) get_par_from_rate_mat(x, index_mat)))
ntips <- do.call(rbind, lapply(full_dat, function(x) length(x$phy$tip.label)))
df_true_long <- get_better_df(df_true, colnames(tmp), "true", ntips)
df_long_list <- list()
for(i in 1:length(df_list)){
  df_long_list[[i]] <- get_better_df(df_list[[i]], colnames(tmp), names(res_list)[i], ntips)
  df_long_list[[i]]$true <- df_true_long$value
  df_long_list[[i]]$diff <- df_long_list[[i]]$value - df_true_long$value
}
df_all <- do.call(rbind, df_long_list)
convert_par_format <- function(par_value) {
  gsub("\\((\\d)\\)_\\((\\d)\\)", "q[\\1][\\2]", par_value)
}
# df_all$par <- sapply(df_all$par, convert_par_format)
df_all_1 <- cbind(sim = "sim-01", df_all)

# sim 2
trees <- lapply(dir("trees/", full.names = TRUE), read.tree)
simulation <- "02"
index_mat <- get_index_mat(nChar=2, nStates=2, nRateClass=1)
# load
par_table_name <- paste0("par_table-", simulation, ".csv")
full_dat_name <- paste0("full_data-", simulation, ".RDS")
cor_dat_name <- paste0("cor_data-", simulation, ".RDS")
sim_res_files <- dir("res/", full.names = TRUE)[grep(paste0("res", simulation, "_"), dir("res/"))]
tmp <- get_par_table(index_mat, 100, mean = 0, sd = 0.25)  
par_table <- read.csv(paste0("parameter_tables/", par_table_name))
colnames(par_table) <- colnames(tmp)
rate_mats <- get_rate_mats(index_mat, par_table)
full_dat <- readRDS(paste0("data/", full_dat_name))
cor_dat <- lapply(full_dat, "[[", "cor_dat")
res_list <- lapply(sim_res_files, readRDS)
names(res_list) <- gsub(".*_", "", sim_res_files) %>% gsub(".RDS", "", .)
# format data and compare results
df_list <- lapply(res_list, function(x) do.call(rbind, lapply(x, get_solution_from_res)))
df_true <- do.call(rbind, lapply(full_dat, function(x) get_par_from_rate_mat(x, index_mat)))
ntips <- do.call(rbind, lapply(full_dat, function(x) length(x$phy$tip.label)))
df_true_long <- get_better_df(df_true, colnames(tmp), "true", ntips)
df_long_list <- list()
for(i in 1:length(df_list)){
  df_long_list[[i]] <- get_better_df(df_list[[i]], colnames(tmp), names(res_list)[i], ntips)
  df_long_list[[i]]$true <- df_true_long$value
  df_long_list[[i]]$diff <- df_long_list[[i]]$value - df_true_long$value
}
df_all <- do.call(rbind, df_long_list)
# convert_par_format <- function(par_value) {
#   gsub("\\((\\d)\\)_\\((\\d)\\)", "q[\\1][\\2]", par_value)
# }
# df_all$par <- sapply(df_all$par, convert_par_format)
df_all_2 <- cbind(sim = "sim-02", df_all)


# sim 3
trees <- lapply(dir("trees/", full.names = TRUE), read.tree)
simulation <- "03"
index_mat <- get_index_mat(nChar=1, nStates=2, nRateClass=2)
# load
par_table_name <- paste0("par_table-", simulation, ".csv")
full_dat_name <- paste0("full_data-", simulation, ".RDS")
cor_dat_name <- paste0("cor_data-", simulation, ".RDS")
sim_res_files <- dir("res/", full.names = TRUE)[grep(paste0("res", simulation, "_"), dir("res/"))]
tmp <- get_par_table(index_mat, 100, mean = 0, sd = 0.25)  
par_table <- read.csv(paste0("parameter_tables/", par_table_name))
colnames(par_table) <- colnames(tmp)
rate_mats <- get_rate_mats(index_mat, par_table)
full_dat <- readRDS(paste0("data/", full_dat_name))
cor_dat <- lapply(full_dat, "[[", "cor_dat")
res_list <- lapply(sim_res_files, readRDS)
names(res_list) <- gsub(".*_", "", sim_res_files) %>% gsub(".RDS", "", .)
# format data and compare results
df_list <- lapply(res_list, function(x) do.call(rbind, lapply(x, get_solution_from_res)))
df_true <- do.call(rbind, lapply(full_dat, function(x) get_par_from_rate_mat(x, index_mat)))
ntips <- do.call(rbind, lapply(full_dat, function(x) length(x$phy$tip.label)))
df_true_long <- get_better_df(df_true, colnames(tmp), "true", ntips)
df_long_list <- list()
for(i in 1:length(df_list)){
  df_long_list[[i]] <- get_better_df(df_list[[i]], colnames(tmp), names(res_list)[i], ntips)
  df_long_list[[i]]$true <- df_true_long$value
  df_long_list[[i]]$diff <- df_long_list[[i]]$value - df_true_long$value
}
df_all <- do.call(rbind, df_long_list)
convert_par_format <- function(par_value) {
  gsub("\\((\\d)\\)_\\((\\d)\\)", "q[\\1][\\2]", par_value)
}
df_all$par <- sapply(df_all$par, convert_par_format)
df_all_3 <- cbind(sim = "sim-02", df_all)
