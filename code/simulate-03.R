# simple two state simulation with hmms
setwd("~/corhmm-dredge/")

library(corHMM)
library(parallel)
library(MASS)

source("code/utils.R")

nSim <- 100
if (detectCores()>100){
  mccores <- 100  
}else{
  mccores <- 4
}

trees <- lapply(dir("trees/", full.names = TRUE), read.tree)
phy <- trees[[1]]

# which simulation number
simulation <- "03"
overwrite <- TRUE

# the various file names
par_table_name <- paste0("par_table-", simulation, ".csv")
full_dat_name <- paste0("full_data-", simulation, ".RDS")
cor_dat_name <- paste0("cor_data-", simulation, ".RDS")

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

# let's create the full data strucutre now.
# a list with phy, sim pars, data
full_dat <- list()
count <- 1
for(i in trees){
  for(j in i){
    for(k in rate_mats){
      tmp <- list(phy = j, par = k, dat = NULL)
      full_dat[[count]] <- tmp
      count <- count+1
    }
  }
}


###### ###### ###### ###### data simulation ###### ###### ###### ###### 
file_found <- full_dat_name %in% dir("data/")
if(!file_found | overwrite){
  for(i in 1:length(full_dat)){
    cat("\r", i, "out of", length(full_dat), "...    ")
    full_dat[[i]]$dat <- get_sim_data(full_dat[[i]]$phy, full_dat[[i]]$par, index_mat)
    full_dat[[i]]$cor_dat <- get_formatted_data(full_dat[[i]]$dat, index_mat)
  }
  saveRDS(full_dat, file = paste0("data/", full_dat_name))
}else{
  full_dat <- readRDS(paste0("data/", full_dat_name))
}

###### ###### ###### ###### model fitting ###### ###### ###### ###### 
file_name <- "res03_unreg.RDS"
file_found <- file_name %in% dir("res/")
if(!file_found | overwrite){
  res_unreg <- mclapply(full_dat, function(x) 
    corHMM(x$phy, x$cor_dat, 2, root.p = "maddfitz"), 
    mc.cores = mccores)
  saveRDS(res_unreg, file = paste0("res/", file_name))
}else{
  res_unreg <- readRDS(paste0("res/", file_name))
}

file_name <- "res03_reg-l1.RDS"
file_found <- file_name %in% dir("res/")
if(!file_found | overwrite){
  res_reg <- mclapply(full_dat, function(x) 
    corHMM:::corHMMDredge(x$phy, x$cor_dat, 2, pen_type = "l1", root.p = "maddfitz"), 
    mc.cores = mccores)
  saveRDS(res_reg, file = paste0("res/", file_name))
}else{
  res_reg <- readRDS(paste0("res/", file_name))
}

file_name <- "res03_reg-l2.RDS"
file_found <- file_name %in% dir("res/")
if(!file_found | overwrite){
  res_reg <- mclapply(full_dat, function(x) 
    corHMM:::corHMMDredge(x$phy, x$cor_dat, 2, pen_type = "l2", root.p = "maddfitz"), 
    mc.cores = mccores)
  saveRDS(res_reg, file = paste0("res/", file_name))
}else{
  res_reg <- readRDS(paste0("res/", file_name))
}

###### ###### ###### ###### summarization ###### ###### ###### ###### 
# supp_hmm <- which((unlist(lapply(res_unreg_b, function(x) x$AICc)) - unlist(lapply(res_unreg, function(x) x$AICc)))>0)
# 
# res_unreg <- res_unreg[unlist(lapply(full_dat, function(x) length(table(x$TipStates)) == 2))]
# res_reg <- res_reg[unlist(lapply(full_dat, function(x) length(table(x$TipStates)) == 2))]
# par_table <- par_table[unlist(lapply(full_dat, function(x) length(table(x$TipStates)) == 2)), ]
# 
# df_unreg <- do.call(rbind, lapply(res_unreg, get_solution_from_res))
# df_reg <- do.call(rbind, lapply(res_reg, get_solution_from_res))
# 
# plot_data <- (cbind(df_unreg, df_reg))
# 
# bias = colMeans(plot_data - cbind(par_table, par_table))
# varr = apply(plot_data, 2, var)
# mse = colMeans((plot_data - cbind(par_table, par_table))^2)
# rmse = sqrt(colMeans((plot_data - cbind(par_table, par_table))^2))
# 
# print(t(data.frame(bias, varr, mse, rmse)))

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
