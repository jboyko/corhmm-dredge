# simple two state simulation with hmms
setwd("~/corhmm-dredge/")

library(corHMM)
library(parallel)
library(MASS)
library(MCMCpack)

source("code/utils.R")

nSim <- 100
mccores <- 100

trees <- lapply(dir("trees/", full.names = TRUE), read.tree)
phy <- trees[[1]]

# which simulation number
simulation <- "03"
overwrite <- FALSE

# the various file names
par_table_name <- paste0("par_table-", simulation, ".csv")
full_dat_name <- paste0("full_data-", simulation, ".RDS")
cor_dat_name <- paste0("cor_data-", simulation, ".RDS")
res_reg_name <- paste0("res_reg-", simulation, ".RDS")
res_unreg_name <- paste0("res_unreg-", simulation, ".RDS")
res_unreg_b_name <- paste0("res_unreg_b-", simulation, ".RDS")
res_bayes_name <- paste0("res_bayes-", simulation, ".RDS")

# creates an index mat appropriate for nchar, nstates, and nhidden
index_mat <- get_index_mat(nChar=1, nStates=2, nRateClass=2)

###### ###### ###### ###### parameter table generation ###### ###### ###### ###### 
file_found <- par_table_name %in% dir("parameter_tables/")
if(!file_found | overwrite){
  par_table <- get_par_table(index_mat, nSim, mean = 0, sd = 0.25)
  # modify for hidden states
  par_table[,c(1,2)] <- par_table[,c(1,2)] * 5
  par_table[,c(3,4)] <- par_table[,c(3,4)] / 5
  par_table[,c(5,6)] <- par_table[,c(5,6)] * 2
  write.csv(par_table, file = paste0("parameter_tables/", par_table_name), row.names = FALSE)
}else{
  tmp <- get_par_table(index_mat, nSim, mean = 0, sd = 0.25)  
  par_table <- read.csv(paste0("parameter_tables/", par_table_name))
  colnames(par_table) <- colnames(tmp)
}

# creates a list of rate matrices for simulation
rate_mats <- get_rate_mats(index_mat, par_table)

###### ###### ###### ###### data simulation ###### ###### ###### ###### 
file_found <- full_dat_name %in% dir("data/")
if(!file_found | overwrite){
  full_dat <- lapply(rate_mats, function(x) get_sim_data(phy, x, index_mat))
  saveRDS(full_dat, file = paste0("data/", full_dat_name))
}else{
  full_dat <- readRDS(paste0("data/", full_dat_name))
}

# format data
cor_dat <- lapply(full_dat, function(x) get_formatted_data(x, index_mat))

###### ###### ###### ###### model fitting ###### ###### ###### ###### 
file_found <- res_unreg_name %in% dir("res/")
if(!file_found | overwrite){
  res_unreg <- mclapply(cor_dat, function(x) 
    corHMM(phy, x, 2), 
    mc.cores = mccores)
  saveRDS(res_unreg, file = paste0("res/", res_unreg_name))
}else{
  res_unreg <- readRDS(paste0("res/", res_unreg_name))
}

file_found <- res_unreg_b_name %in% dir("res/")
if(!file_found | overwrite){
  res_unreg_b <- mclapply(cor_dat, function(x) 
    corHMM(phy, x, 1), 
    mc.cores = mccores)
  saveRDS(res_unreg_b, file = paste0("res/", res_unreg_b_name))
}else{
  res_unreg_b <- readRDS(paste0("res/", res_unreg_b_name))
}


file_found <- res_reg_name %in% dir("res/")
if(!file_found | overwrite){
  res_reg <- mclapply(cor_dat, function(x) 
    corHMM:::corHMMDredge(phy, x, 2, pen_type = "logl1"), 
    mc.cores = mccores)
  saveRDS(res_reg, file = paste0("res/", res_reg_name))
}else{
  res_reg <- readRDS(paste0("res/", res_reg_name))
}

file_found <- res_bayes_name %in% dir("res/")
if(!file_found | overwrite){
  nPar <- max(index_mat$full_rate_mat)
  res_bayes <- mclapply(cor_dat, function(x) 
    MCMCmetrop1R(log_posterior, theta.init=rep(0.5, nPar), force.samp=TRUE,
                 optim.lower = 1e-8, optim.upper = 1e2, optim.method = "L-BFGS-B",
                 mcmc=10000, burnin=1000, thin=10, verbose=TRUE, 
                 tree=phy, data=x, rate.cat = 2, logfun=TRUE),
    mc.cores = mccores)
  saveRDS(res_bayes, file = paste0("res/", res_bayes_name))
}else{
  res_bayes <- readRDS(paste0("res/", res_bayes_name))
}

###### ###### ###### ###### summarization ###### ###### ###### ###### 
supp_hmm <- which((unlist(lapply(res_unreg_b, function(x) x$AICc)) - unlist(lapply(res_unreg, function(x) x$AICc)))>0)

res_unreg <- res_unreg[unlist(lapply(full_dat, function(x) length(table(x$TipStates)) == 2))]
res_reg <- res_reg[unlist(lapply(full_dat, function(x) length(table(x$TipStates)) == 2))]
par_table <- par_table[unlist(lapply(full_dat, function(x) length(table(x$TipStates)) == 2)), ]

df_unreg <- do.call(rbind, lapply(res_unreg, get_solution_from_res))
df_reg <- do.call(rbind, lapply(res_reg, get_solution_from_res))

plot_data <- (cbind(df_unreg, df_reg))

bias = colMeans(plot_data - cbind(par_table, par_table))
varr = apply(plot_data, 2, var)
mse = colMeans((plot_data - cbind(par_table, par_table))^2)
rmse = sqrt(colMeans((plot_data - cbind(par_table, par_table))^2))

print(t(data.frame(bias, varr, mse, rmse)))

# boxplot(plot_data); abline(h = colMeans(par_table))


# bias = colMeans(log(plot_data) - log(cbind(par_table, par_table)))
# varr = apply(log(plot_data), 2, var)
# mse = colMeans((log(plot_data) - log(cbind(par_table, par_table)))^2)
# rmse = sqrt(colMeans((log(plot_data) - log(cbind(par_table, par_table)))^2))

# t(data.frame(bias, varr, mse, rmse))

# boxplot(log(plot_data)); abline(h = colMeans(log(par_table)))

# colMeans(par_table)
# colMeans(plot_data)

# df_reg[supp_hmm,]
