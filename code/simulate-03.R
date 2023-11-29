# simple two state simulation with hmms
setwd("~/corhmm-dredge/")

library(corHMM)
library(parallel)
library(MASS)
source("code/utils.R")

# creates an index mat appropriate for nchar, nstates, and nhidden
index_mat <- get_index_mat(nChar=1, nStates=2, nRateClass=2)

# samples the possible parameter values 
par_table <- get_par_table(index_mat, nSim, mean = 0, sd = 0.25)

# modify for hidden states
par_table[,c(1,2)] <- par_table[,c(1,2)] * 5
par_table[,c(3,4)] <- par_table[,c(3,4)] / 5
par_table[,c(5,6)] <- par_table[,c(5,6)] * 2

# creates a list of rate matrices for simulation
rate_mats <- get_rate_mats(index_mat, par_table)

# simulate data
full_dat <- lapply(rate_mats, function(x) get_sim_data(trees[[1]], x, index_mat))

# format data
cor_dat <- lapply(full_dat, function(x) get_formatted_data(x, index_mat))

# fit model to data
res_unreg <- mclapply(cor_dat, function(x) 
  corHMM(phy, x, 2), 
  mc.cores = 10)

res_unreg_b <- mclapply(cor_dat, function(x) 
  corHMM(phy, x, 1), 
  mc.cores = 10)

supp_hmm <- which((unlist(lapply(res_unreg_b, function(x) x$AICc)) - unlist(lapply(res_unreg, function(x) x$AICc)))>-2)

res_reg <- mclapply(cor_dat, function(x) 
  corHMM:::corHMMDredge(phy, x, 2, pen_type = "logl1", lambda = 1), 
  mc.cores = 10)


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
