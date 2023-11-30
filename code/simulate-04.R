# simple two state simulation with hmms
setwd("~/corhmm-dredge/")

library(corHMM)
library(parallel)
library(MASS)
source("code/utils.R")

nSim = 100
trees <- lapply(dir("trees/", full.names = TRUE), read.tree)
phy <- trees[[1]]

# creates an index mat appropriate for nchar, nstates, and nhidden
index_mat <- get_index_mat(nChar=4, nStates=2, nRateClass=1)

# samples the possible parameter values 
par_table <- get_par_table(index_mat, nSim, mean = 0, sd = 0.25)

# creates a list of rate matrices for simulation
rate_mats <- get_rate_mats(index_mat, par_table)

# simulate data
full_dat <- lapply(rate_mats, function(x) get_sim_data(trees[[1]], x, index_mat))

# lapply(full_dat, function(x) table(x$TipStates))

# format data
cor_dat <- lapply(full_dat, function(x) get_formatted_data(x, index_mat))

# fit model to data
res_unreg <- mclapply(cor_dat, function(x) 
  corHMM(phy, x, 1), 
  mc.cores = 10)

res_reg <- mclapply(cor_dat, function(x) 
  corHMM:::corHMMDredge(phy, x, 1, pen_type = "logl1", lambda = 1), 
  mc.cores = 10)

res_unreg <- res_unreg[unlist(lapply(full_dat, function(x) length(table(x$TipStates)) == 8))]
res_reg <- res_reg[unlist(lapply(full_dat, function(x) length(table(x$TipStates)) == 8))]
par_table <- par_table[unlist(lapply(full_dat, function(x) length(table(x$TipStates)) == 8)), ]

df_unreg <- do.call(rbind, lapply(res_unreg, get_solution_from_res))
df_reg <- do.call(rbind, lapply(res_reg, get_solution_from_res))

plot_data <- (cbind(df_unreg, df_reg))

bias = colMeans(plot_data - cbind(par_table, par_table))
varr = apply(plot_data, 2, var)
mse = colMeans((plot_data - cbind(par_table, par_table))^2)
rmse = sqrt(colMeans((plot_data - cbind(par_table, par_table))^2))

t(data.frame(bias, varr, mse, rmse))

boxplot(plot_data); abline(h = colMeans(par_table))


bias = colMeans(log(plot_data) - log(cbind(par_table, par_table)))
varr = apply(log(plot_data), 2, var)
mse = colMeans((log(plot_data) - log(cbind(par_table, par_table)))^2)
rmse = sqrt(colMeans((log(plot_data) - log(cbind(par_table, par_table)))^2))

t(data.frame(bias, varr, mse, rmse))

boxplot(log(plot_data)); abline(h = colMeans(log(par_table)))

colMeans(par_table)
colMeans(plot_data)

df_reg[supp_hmm,]
