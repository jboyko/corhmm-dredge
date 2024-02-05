setwd("~/corhmm-dredge/")

library(corHMM)
library(parallel)
library(MASS)
library(dplyr)
library(ggplot2)
library(tidyr)

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
df_unreg <- do.call(rbind, lapply(res_list[[2]], get_solution_from_res))
df_list <- lapply(res_list, function(x) do.call(rbind, lapply(x, get_solution_from_res)))
df_unreg <- do.call(rbind, lapply(res_unreg, get_solution_from_res))
df_reg <- do.call(rbind, lapply(res_reg, get_solution_from_res))
df_true <- do.call(rbind, lapply(full_dat, function(x) get_par_from_rate_mat(x, index_mat)))
ntips <- do.call(rbind, lapply(full_dat, function(x) length(x$phy$tip.label)))
# df_bayes <- do.call(rbind, lapply(res_bayes, get_solution_from_res))

df_true_long <- get_better_df(df_true, colnames(tmp), "true", ntips)
df_unreg_long <- get_better_df(df_unreg, colnames(tmp), "unreg", ntips)
df_reg_long <- get_better_df(df_reg, colnames(tmp), "reg", ntips)

df_unreg_long$diff <- df_unreg_long$value - df_true_long$value
df_reg_long$diff <- df_reg_long$value - df_true_long$value
  
# Calculate the differences between df_unreg and df_true
# df_unreg_diff <- df_unreg %>%
#   left_join(df_true, by = c("ntips", "par")) %>%
#   mutate(diff_unreg = value.x - value.y)

# Calculate MSE and RMSE for df_reg_diff
df_reg_summary <- df_reg_long %>%
  group_by(par, ntips) %>%
  summarize(
    mean_diff_reg = mean(diff),
    var_diff_reg = var(diff),
    mse_reg = mean(diff^2),  # Calculate MSE
    rmse_reg = sqrt(mse_reg)     # Calculate RMSE
  )

# Calculate MSE and RMSE for df_unreg_diff
df_unreg_summary <- df_unreg_long %>%
  group_by(par, ntips) %>%
  summarize(
    mean_diff_unreg = mean(diff),
    var_diff_unreg = var(diff),
    mse_unreg = mean(diff^2),  # Calculate MSE
    rmse_unreg = sqrt(mse_unreg)     # Calculate RMSE
  )

# Merge the summaries for df_reg and df_unreg
comparison_summary <- df_reg_summary %>%
  left_join(df_unreg_summary, by = c("par", "ntips"))

print(comparison_summary)



ggplot(df_reg_long, aes(x = value, fill = par)) +
  geom_histogram(binwidth = 0.1, position = "identity", alpha = 0.5) +
  facet_wrap(~par) +
  labs(title = "Comparison of df_reg_long and df_unreg_long")

ggplot(df_unreg_long, aes(x = value, fill = par)) +
  geom_histogram(binwidth = 0.1, position = "identity", alpha = 0.5) +
  facet_wrap(~par) +
  labs(title = "Comparison of df_unreg_long and df_reg_long")

# Perform t-test for each level of par
t_test_results <- df_reg_long %>%
  group_by(par) %>%
  summarize(p_value = t.test(value, df_unreg_long$value)$p.value)



par(mfrow=c(2,2))
plot(log(df_true[,1]), log(df_unreg[,1]))
abline(coef = c(0,1), col = "red")
plot(df_true[,2], df_unreg[,2])
plot(log(df_true[,1]), log(df_reg[,1]))
abline(coef = c(0,1), col = "red")
plot(df_true[,2], df_reg[,2])
dev.off()

# plot_data <- (cbind(df_unreg, df_reg, df_bayes))
# colnames(plot_data) <- paste0(colnames(tmp), rep(c("-unreg", "-reg", "-bayes"), each = 2))

plot_data <- (cbind(df_unreg, df_reg))
colnames(plot_data) <- paste0(colnames(tmp), rep(c("-unreg", "-reg"), each = 2))

bias = colMeans(plot_data - cbind(df_true, df_true), na.rm = TRUE)
varr = apply(plot_data - cbind(df_true, df_true), 2, function(x) var(x, na.rm = TRUE))
mse = colMeans((plot_data - cbind(df_true, df_true))^2, na.rm = TRUE)
rmse = sqrt(colMeans((plot_data - cbind(df_true,df_true))^2, na.rm = TRUE))

print((data.frame(bias, varr, mse, rmse)))

plot_data_long <- pivot_longer(as.data.frame(plot_data), cols = everything())
plot_data_long <- data.frame(do.call(rbind, strsplit(plot_data_long$name, "-")), value = plot_data_long$value)
colnames(plot_data_long) <- c("trans", "type", "value")

ggplot(data = plot_data_long, aes(x = type, y = log(value))) +
  geom_boxplot() +
  theme_minimal() +
  facet_wrap(~trans)

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

