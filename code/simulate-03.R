# simple 3 character simulation with binary states
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
index_mat <- get_index_mat(nChar=3, nStates=2, nRateClass=1)

###### ###### ###### ###### parameter table generation ###### ###### ###### ###### 
file_found <- par_table_name %in% dir("parameter_tables/")
if(!file_found | overwrite){
  par_table <- get_par_table(index_mat, nSim, mean = 0, sd = 0.25)  
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
    tip_states <- rownames(index_mat$full_rate_mat)[full_dat[[i]]$dat$TipStates]
    tmp_dat <- do.call(rbind, strsplit(tip_states, "|"))[,c(1,3,5)]
    cor_dat <- data.frame(sp = names(full_dat[[i]]$dat$TipStates),
                          x1 = factor(tmp_dat[,1], c(1,2)),
                          x2 = factor(tmp_dat[,2], c(1,2)),
                          x3 = factor(tmp_dat[,3], c(1,2)))
    full_dat[[i]]$cor_dat <- cor_dat
  }
  saveRDS(full_dat, file = paste0("data/", full_dat_name))
}else{
  full_dat <- readRDS(paste0("data/", full_dat_name))
}

###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ######
###### ###### ###### ###### model fitting ###### ###### ###### ###### ###### ######
###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ######

results_dir <- get_full_path("param_results/")
fit_models(full_dat, "l1", 0, TRUE, paste0(results_dir, "/res03_l0.RDS"), overwrite, mccores)
fit_models(full_dat, "l1", 1, TRUE, paste0(results_dir, "/res03_l1.RDS"), overwrite, mccores)
fit_models(full_dat, "l2", 1, TRUE, paste0(results_dir, "/res03_l2.RDS"), overwrite, mccores)
fit_models(full_dat, "er", 1, TRUE, paste0(results_dir, "/res03_er.RDS"), overwrite, mccores)
