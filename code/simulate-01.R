# simple binary simulation
setwd("~/corhmm-dredge/")

library(corHMM)
library(parallel)
library(MASS)
source("code/utils.R")

nSim = 100
trees <- lapply(dir("trees/", full.names = TRUE), read.tree)
phy <- trees[[1]]

# creates an index mat appropriate for nchar, nstates, and nhidden
index_mat <- get_index_mat(nChar=1, nStates=2, nRateClass=1)

# samples the possible parameter values 
par_table <- get_par_table(index_mat, nSim, mean = 0, sd = 0.25)

# creates a list of rate matrices for simulation
rate_mats <- get_rate_mats(index_mat, par_table)
rate_mat = rate_mats[[1]]
# simulate data
full_dat <- lapply(rate_mats, function(x) get_sim_data(trees[[1]], x, index_mat))

# format data
cor_dat <- lapply(full_dat, function(x) get_formatted_data(x, index_mat))

# fit model to data
res_unreg <- mclapply(cor_dat, function(x) 
  corHMM(phy, x, 1), 
  mc.cores = 10)

res_reg <- mclapply(cor_dat, function(x) 
  corHMM:::corHMMDredge(phy, x, 1, pen_type = "logl1"), 
  mc.cores = 10)


df_unreg <- do.call(rbind, lapply(res_unreg, get_solution_from_res))
df_reg <- do.call(rbind, lapply(res_reg, get_solution_from_res))

plot_data <- (cbind(df_unreg, df_reg))
colnames(plot_data) <- c("01_unreg", "10_unreg", "01_reg", "10_reg")

bias = colMeans(plot_data) - 1
varr = apply(plot_data, 2, var)
mse = colMeans((plot_data - 1)^2)
rmse = sqrt(colMeans((plot_data - 1)^2))

t(data.frame(bias, varr, mse, rmse))

boxplot(plot_data); abline(h = 1)

