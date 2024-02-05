setwd("~/corhmm-dredge/")

library(corHMM)
library(parallel)
library(MASS)
library(dplyr)
library(ggplot2)
library(tidyr)
library(gghalves)

source("code/utils.R")

trees <- lapply(dir("trees/", full.names = TRUE), read.tree)
phy <- trees[[1]]
nSim = 100

############### Simulation scenario 1 ####################
# which simulation number
simulation <- "01"

# the various file names
par_table_name <- paste0("par_table-", simulation, ".csv")
full_dat_name <- paste0("full_data-", simulation, ".RDS")
cor_dat_name <- paste0("cor_data-", simulation, ".RDS")

sim_res_files <- dir("res/", full.names = TRUE)[grep(paste0("res", simulation, "_"), dir("res/"))]

# res_reg_name <- paste0("res_reg-", simulation, ".RDS")
# res_unreg_name <- paste0("res_unreg-", simulation, ".RDS")
# res_bayes_name <- paste0("res_bayes-", simulation, ".RDS")

# load everything
index_mat <- get_index_mat(nChar=1, nStates=2, nRateClass=1)
tmp <- get_par_table(index_mat, nSim, mean = 0, sd = 0.25)  
par_table <- read.csv(paste0("parameter_tables/", par_table_name))
colnames(par_table) <- colnames(tmp)
rate_mats <- get_rate_mats(index_mat, par_table)
full_dat <- readRDS(paste0("data/", full_dat_name))
cor_dat <- lapply(full_dat, "[[", "cor_dat")

# load results
res_list <- lapply(sim_res_files, readRDS)
names(res_list) <- gsub(".*_", "", sim_res_files) %>% gsub(".RDS", "", .)

# format data and compare results
df_list <- lapply(res_list, function(x) do.call(rbind, lapply(x, get_solution_from_res)))
df_true <- do.call(rbind, lapply(full_dat, function(x) get_par_from_rate_mat(x, index_mat)))
ntips <- do.call(rbind, lapply(full_dat, function(x) length(x$phy$tip.label)))
# df_bayes <- do.call(rbind, lapply(res_bayes, get_solution_from_res))

df_true_long <- get_better_df(df_true, colnames(tmp), "true", ntips)
df_long_list <- list()
for(i in 1:length(df_list)){
  df_long_list[[i]] <- get_better_df(df_list[[i]], colnames(tmp), names(res_list)[i], ntips)
  df_long_list[[i]]$diff <- df_long_list[[i]]$value - df_true_long$value
}

df_all <- do.call(rbind, df_long_list)

# Calculate MSE and RMSE for df_reg_diff
df_summary <- df_all %>%
  group_by(ntips, type, par) %>%
  summarize(
    bias = mean(diff),
    var = var(diff),
    mse= mean(diff^2),  # Calculate MSE
    rmse = sqrt(mse)     # Calculate RMSE
  )

print(df_summary)

ggplot(df_all, aes(x = factor(par), y = diff, fill = type, color = type)) +
  geom_violin() +
  facet_grid(ntips~.) +
  coord_cartesian(ylim=c(-5, 5))


n=0.02
ggplot(df_all, aes(x = factor(type), y = diff, color = type)) + 
  geom_half_boxplot(center=TRUE, errorbar.draw=FALSE, width=0.75, nudge=n, outlier.colour = NA) +
  geom_half_point(alpha = 0.25) +
  geom_half_violin(side="r", nudge=n) +
  coord_cartesian(ylim=c(-5,5)) +
  facet_wrap(par~ntips) +
  theme_bw()

ggsave("plots/plot_01.pdf")

# ASR COMPARISON
max_states_unreg <- lapply(res_unreg, function(x) apply(x$states, 1, which.max))
max_states_reg <- lapply(res_reg, function(x) apply(x$states, 1, which.max))
true_states <- lapply(full_dat, function(x) x$dat$NodeStates)

table(mapply(function(x, y) x == y, x = max_states_unreg, y = true_states))
table(mapply(function(x, y) x == y, x = max_states_reg, y = true_states))

############### Simulation scenario 2 ####################
# which simulation number
simulation <- "02"

# the various file names
par_table_name <- paste0("par_table-", simulation, ".csv")
full_dat_name <- paste0("full_data-", simulation, ".RDS")
cor_dat_name <- paste0("cor_data-", simulation, ".RDS")
res_reg_name <- paste0("res_reg-", simulation, ".RDS")
res_unreg_name <- paste0("res_unreg-", simulation, ".RDS")
# res_bayes_name <- paste0("res_bayes-", simulation, ".RDS")

# load everything
index_mat <- get_index_mat(nChar=2, nStates=2, nRateClass=1)
tmp <- get_par_table(index_mat, nSim, mean = 0, sd = 0.25)  
par_table <- read.csv(paste0("parameter_tables/", par_table_name))
colnames(par_table) <- colnames(tmp)
rate_mats <- get_rate_mats(index_mat, par_table)
full_dat <- readRDS(paste0("data/", full_dat_name))
cor_dat <- lapply(full_dat, function(x) get_formatted_data(x, index_mat))
res_unreg <- readRDS(paste0("res/", res_unreg_name))
res_reg <- readRDS(paste0("res/", res_reg_name))
# res_bayes <- readRDS(paste0("res/", res_bayes_name))

df_unreg <- do.call(rbind, lapply(res_unreg, get_solution_from_res))
df_reg <- do.call(rbind, lapply(res_reg, get_solution_from_res))
# df_bayes <- do.call(rbind, lapply(res_bayes, get_solution_from_res))

max_states_unreg <- lapply(res_unreg, function(x) apply(x$states, 1, which.max))
max_states_reg <- lapply(res_reg, function(x) apply(x$states, 1, which.max))
true_states <- lapply(full_dat, "[[", "NodeStates")

table(mapply(function(x, y) x == y, x = max_states_unreg, y = true_states))
table(mapply(function(x, y) x == y, x = max_states_reg, y = true_states))

# plot_data <- (cbind(df_unreg, df_reg, df_bayes))
# colnames(plot_data) <- paste0(colnames(tmp), rep(c("-unreg", "-reg", "-bayes"), each = length(colnames(tmp))))

plot_data <- (cbind(df_unreg, df_reg))
colnames(plot_data) <- paste0(colnames(tmp), rep(c("-unreg", "-reg"), each = length(colnames(tmp))))

bias = colMeans(plot_data - cbind(par_table, par_table), na.rm = TRUE)
varr = apply(plot_data, 2, function(x) var(x, na.rm = TRUE))
mse = colMeans((plot_data - cbind(par_table, par_table))^2, na.rm = TRUE)
rmse = sqrt(colMeans((plot_data - cbind(par_table,par_table))^2, na.rm = TRUE))

print(t(data.frame(bias, varr, mse, rmse)))

plot_data_long <- pivot_longer(as.data.frame(plot_data), cols = everything())
plot_data_long <- data.frame(do.call(rbind, strsplit(plot_data_long$name, "-")), value = plot_data_long$value)
colnames(plot_data_long) <- c("trans", "type", "value")

ggplot(data = plot_data_long, aes(x = type, y = log(value))) +
  geom_boxplot() +
  theme_minimal() +
  facet_wrap(~trans)

############### Simulation scenario 3 ####################
# which simulation number
simulation <- "03"

# the various file names
par_table_name <- paste0("par_table-", simulation, ".csv")
full_dat_name <- paste0("full_data-", simulation, ".RDS")
cor_dat_name <- paste0("cor_data-", simulation, ".RDS")
res_reg_name <- paste0("res_reg-", simulation, ".RDS")
res_unreg_name <- paste0("res_unreg-", simulation, ".RDS")
# res_bayes_name <- paste0("res_bayes-", simulation, ".RDS")

# load everything
index_mat <- get_index_mat(nChar=1, nStates=2, nRateClass=2)
tmp <- get_par_table(index_mat, nSim, mean = 0, sd = 0.25)  
par_table <- read.csv(paste0("parameter_tables/", par_table_name))
colnames(par_table) <- colnames(tmp)
rate_mats <- get_rate_mats(index_mat, par_table)
full_dat <- readRDS(paste0("data/", full_dat_name))
cor_dat <- lapply(full_dat, function(x) get_formatted_data(x, index_mat))
res_unreg <- readRDS(paste0("res/", res_unreg_name))
res_reg <- readRDS(paste0("res/", res_reg_name))
# res_bayes <- readRDS(paste0("res/", res_bayes_name))

df_unreg <- do.call(rbind, lapply(res_unreg, get_solution_from_res))
df_reg <- do.call(rbind, lapply(res_reg, get_solution_from_res))
# df_bayes <- do.call(rbind, lapply(res_bayes, get_solution_from_res))

# plot_data <- (cbind(df_unreg, df_reg, df_bayes))
# colnames(plot_data) <- paste0(colnames(tmp), rep(c("-unreg", "-reg", "-bayes"), each = length(colnames(tmp))))

plot_data <- (cbind(df_unreg, df_reg))
colnames(plot_data) <- paste0(colnames(tmp), rep(c("-unreg", "-reg"), each = length(colnames(tmp))))

bias = colMeans(plot_data - cbind(par_table, par_table), na.rm = TRUE)
varr = apply(plot_data, 2, function(x) var(x, na.rm = TRUE))
mse = colMeans((plot_data - cbind(par_table, par_table))^2, na.rm = TRUE)
rmse = sqrt(colMeans((plot_data - cbind(par_table,par_table))^2, na.rm = TRUE))

print(t(data.frame(bias, varr, mse, rmse)))

plot_data_long <- pivot_longer(as.data.frame(plot_data), cols = everything())
plot_data_long <- data.frame(do.call(rbind, strsplit(plot_data_long$name, "-")), value = plot_data_long$value)
colnames(plot_data_long) <- c("trans", "type", "value")

ggplot(data = plot_data_long, aes(x = type, y = log(value))) +
  geom_boxplot() +
  theme_minimal() +
  facet_wrap(~trans)

max_states_unreg <- lapply(res_unreg, function(x) apply(x$states, 1, which.max))
max_states_reg <- lapply(res_reg, function(x) apply(x$states, 1, which.max))
true_states <- lapply(full_dat, "[[", "NodeStates")

table(mapply(function(x, y) x == y, x = max_states_unreg, y = true_states))
table(mapply(function(x, y) x == y, x = max_states_reg, y = true_states))


############### Simulation scenario 4 ####################
# which simulation number
simulation <- "04"

# the various file names
par_table_name <- paste0("par_table-", simulation, ".csv")
full_dat_name <- paste0("full_data-", simulation, ".RDS")
cor_dat_name <- paste0("cor_data-", simulation, ".RDS")
res_reg_name <- paste0("res_reg-", simulation, ".RDS")
res_unreg_name <- paste0("res_unreg-", simulation, ".RDS")
# res_bayes_name <- paste0("res_bayes-", simulation, ".RDS")

# load everything
index_mat <- get_index_mat(nChar=3, nStates=2, nRateClass=1)
tmp <- get_par_table(index_mat, nSim, mean = 0, sd = 0.25)  
par_table <- read.csv(paste0("parameter_tables/", par_table_name))
colnames(par_table) <- colnames(tmp)
rate_mats <- get_rate_mats(index_mat, par_table)
full_dat <- readRDS(paste0("data/", full_dat_name))
cor_dat <- lapply(full_dat, function(x) get_formatted_data(x, index_mat))
res_unreg <- readRDS(paste0("res/", res_unreg_name))
res_reg <- readRDS(paste0("res/", res_reg_name))
# res_bayes <- readRDS(paste0("res/", res_bayes_name))

df_unreg <- do.call(rbind, lapply(res_unreg, function(x) 
  get_solution_from_res(x, index_mat$full_rate_mat)))
df_reg <- do.call(rbind, lapply(res_reg, function(x) 
  get_solution_from_res(x, index_mat$full_rate_mat)))
# df_bayes <- do.call(rbind, lapply(res_bayes, get_solution_from_res))

# plot_data <- (cbind(df_unreg, df_reg, df_bayes))
# colnames(plot_data) <- paste0(colnames(tmp), rep(c("-unreg", "-reg", "-bayes"), each = length(colnames(tmp))))

plot_data <- (cbind(df_unreg, df_reg))
colnames(plot_data) <- paste0(colnames(tmp), rep(c("-unreg", "-reg"), each = length(colnames(tmp))))

bias = colMeans(plot_data - cbind(par_table, par_table), na.rm = TRUE)
varr = apply(plot_data, 2, function(x) var(x, na.rm = TRUE))
mse = colMeans((plot_data - cbind(par_table, par_table))^2, na.rm = TRUE)
rmse = sqrt(colMeans((plot_data - cbind(par_table,par_table))^2, na.rm = TRUE))

print(t(data.frame(bias, varr, mse, rmse)))

plot_data_long <- pivot_longer(as.data.frame(plot_data), cols = everything())
plot_data_long <- data.frame(do.call(rbind, strsplit(plot_data_long$name, "-")), value = plot_data_long$value)
colnames(plot_data_long) <- c("trans", "type", "value")

ggplot(data = plot_data_long, aes(x = type, y = log(value))) +
  geom_violin() +
  theme_minimal() +
  facet_wrap(~trans)

max_states_unreg <- lapply(res_unreg, function(x) apply(x$states, 1, which.max))
max_states_reg <- lapply(res_reg, function(x) apply(x$states, 1, which.max))
true_states <- lapply(full_dat, "[[", "NodeStates")

table(mapply(function(x, y) x == y, x = max_states_unreg, y = true_states))
table(mapply(function(x, y) x == y, x = max_states_reg, y = true_states))

