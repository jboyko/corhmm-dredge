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
res_reg_name <- paste0("res_reg-", simulation, ".RDS")
res_unreg_name <- paste0("res_unreg-", simulation, ".RDS")
res_bayes_name <- paste0("res_bayes-", simulation, ".RDS")

# load everything
index_mat <- get_index_mat(nChar=1, nStates=2, nRateClass=1)
tmp <- get_par_table(index_mat, nSim, mean = 0, sd = 0.25)  
par_table <- read.csv(paste0("parameter_tables/", par_table_name))
colnames(par_table) <- colnames(tmp)
rate_mats <- get_rate_mats(index_mat, par_table)
full_dat <- readRDS(paste0("data/", full_dat_name))
cor_dat <- lapply(full_dat, function(x) get_formatted_data(x, index_mat))
res_unreg <- readRDS(paste0("res/", res_unreg_name))
res_reg <- readRDS(paste0("res/", res_reg_name))

df_unreg <- do.call(rbind, lapply(res_unreg, get_solution_from_res))
df_reg <- do.call(rbind, lapply(res_reg, get_solution_from_res))

plot_data <- (cbind(df_unreg, df_reg))
colnames(plot_data) <- c("01_unreg", "10_unreg", "01_reg", "10_reg")

bias = colMeans(plot_data - cbind(par_table, par_table))
varr = apply(plot_data, 2, var)
mse = colMeans((plot_data - cbind(par_table, par_table))^2)
rmse = sqrt(colMeans((plot_data - cbind(par_table,par_table))^2))

print(t(data.frame(bias, varr, mse, rmse)))

plot_data_long <- pivot_longer(as.data.frame(plot_data), cols = everything())
plot_data_long <- data.frame(do.call(rbind, strsplit(plot_data_long$name, "_")), value = plot_data_long$value)
colnames(plot_data_long) <- c("trans", "type", "value")

ggplot(data = plot_data_long, aes(x = type, y = value))

############### Simulation scenario 2 ####################

