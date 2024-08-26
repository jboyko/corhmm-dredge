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
df_out <- list()
count <- 0
for (simulation in c("01", "02", "04")){
  print(simulation)
  count <- count + 1
  # the various file names
  par_table_name <- paste0("par_table-", simulation, ".csv")
  full_dat_name <- paste0("full_data-", simulation, ".RDS")
  cor_dat_name <- paste0("cor_data-", simulation, ".RDS")
  
  sim_res_files <- dir("results/", full.names = TRUE)[grep(paste0("res", simulation), dir("results/"))]
  
  # load everything
  if(simulation == "01"){
    index_mat <- get_index_mat(nChar=1, nStates=2, nRateClass=1)
  }
  if(simulation == "02"){
    index_mat <- get_index_mat(nChar=2, nStates=2, nRateClass=1)
  }
  if(simulation == "03"){
    index_mat <- get_index_mat(nChar=1, nStates=2, nRateClass=2)
  }
  if(simulation == "04"){
    index_mat <- get_index_mat(nChar=3, nStates=2, nRateClass=1)
  }
  
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
  df_list <- lapply(res_list, function(x) 
    do.call(rbind, lapply(x, function(y) 
      get_solution_from_res(y, index_mat$rate_mat))))
  # index_na <- do.call(cbind, lapply(df_list, function(y) apply(y, 1, function(x) any(is.na(x)))))
  # for(i in 1:length(df_list)){
  #   df_list[[i]] <- df_list[[i]][!apply(index_na, 1, any),]
  # }
  df_true <- do.call(rbind, lapply(full_dat, function(x) get_par_from_rate_mat(x, index_mat)))
  ntips <- do.call(rbind, lapply(full_dat, function(x) length(x$phy$tip.label)))
  # df_bayes <- do.call(rbind, lapply(res_bayes, get_solution_from_res))
  
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
  
  # Calculate MSE and RMSE for df_reg_diff
  df_summary <- df_all %>%
    group_by(ntips, type, par) %>%
    summarize(
      bias = mean(diff, na.rm=TRUE),
      var = var(diff, na.rm=TRUE),
      mse= mean(diff^2, na.rm=TRUE),  # Calculate MSE
      rmse = sqrt(mse)     # Calculate RMSE
    )
  df_out[[count]] <- df_all
}

ggplot(df_out[[1]], 
  aes(x = factor(type), y = diff, fill = as.factor(ntips))) + 
  geom_boxplot(width=0.75, outlier.colour = NA) +
  coord_cartesian(ylim=c(-5, 5))

# Define the list of functions
agg_funcs <- list(mean = mean, varriance = var, median = median)

# with 100 trans / MY
dat <- do.call(rbind, df_out)
dat$rmse <-(dat$value - dat$true)^2
# Apply the aggregate function
result <- aggregate(dat[,c(4:7)], 
  by = list(ntips = dat$ntips, type = dat$type), 
  FUN = function(x) sapply(agg_funcs, function(f) f(x, na.rm = TRUE))
)
result <- cbind(include="all",result)
# Unlist and convert to a data frame for better readability
final_df <- data.frame(result)

# <50 trans / MY
dat_prune_50 <- as.data.frame(dat)[!dat$value > 50,]
dat_prune_50$rmse <-(dat_prune_50$value - dat_prune_50$true)^2
# Apply the aggregate function
result <- aggregate(dat_prune_50[,4:7], 
  by = list(ntips = dat_prune_50$ntips, type = dat_prune_50$type), 
  FUN = function(x) sapply(agg_funcs, function(f) f(x, na.rm = TRUE))
)
result <- cbind(include="<50",result)
result$value[7:9,2]/result$value[,2]
final_df <- rbind(final_df, result)

# <10 trans / MY
dat_prune_10 <- as.data.frame(dat)[!dat$value > 10,]
dat_prune_10$rmse <-(dat_prune_10$value - dat_prune_10$true)^2
# Apply the aggregate function
result <- aggregate(dat_prune_10[,4:7], 
  by = list(ntips = dat_prune_10$ntips, type = dat_prune_10$type), 
  FUN = function(x) sapply(agg_funcs, function(f) f(x, na.rm = TRUE))
)
result <- cbind(include="<10",result)
result$value[7:9,2]/result$value[,2]
final_df <- rbind(final_df, result)
final_df$rmse <- sqrt(final_df$rmse)

colnames(final_df) <- c("include", "ntips", "type", "estimate", "sim", "diff", "rmse")

write.csv(final_df, "tables/par-summary-by-type.csv", row.names = FALSE)

stats_df <- data.frame(
  include = final_df$include,
  type = final_df$type,
  ntips = final_df$ntips,
  bias=final_df$diff[,1],
  varriance=final_df$estimate[,2],
  rmse=final_df$rmse[,1])

write.csv(stats_df, "tables/stats-summary-by-type.csv", row.names = FALSE)


# with 100 trans / MY
dat <- do.call(rbind, df_out)
# Apply the aggregate function
result <- aggregate(dat[,4:6], 
  by = list(ntips = dat$ntips, 
    type = dat$type, 
    par = dat$par), 
  FUN = function(x) sapply(agg_funcs, function(f) f(x, na.rm = TRUE))
)
result <- cbind(include="all",result)
# Unlist and convert to a data frame for better readability
final_df <- data.frame(result)

# <50 trans / MY
dat_prune_50 <- as.data.frame(dat)[!dat$value > 50,]
# Apply the aggregate function
result <- aggregate(dat_prune_50[,4:6], 
  by = list(ntips = dat_prune_50$ntips, 
    type = dat_prune_50$type, 
    par = dat_prune_50$par), 
  FUN = function(x) sapply(agg_funcs, function(f) f(x, na.rm = TRUE))
)
result <- cbind(include="<50",result)
result$value[7:9,2]/result$value[,2]
final_df <- rbind(final_df, result)

# <10 trans / MY
dat_prune_10 <- as.data.frame(dat)[!dat$value > 10,]
# Apply the aggregate function
result <- aggregate(dat_prune_10[,4:6], 
  by = list(ntips = dat_prune_10$ntips, 
    type = dat_prune_10$type,
    par = dat_prune_10$par), 
  FUN = function(x) sapply(agg_funcs, function(f) f(x, na.rm = TRUE))
)
result <- cbind(include="<10",result)
result$value[7:9,2]/result$value[,2]
final_df <- rbind(final_df, result)

colnames(final_df) <- c("include", "ntips", "type", "par", "estimate", "sim", "diff")

write.csv(final_df, "tables/par-summary-by-type-par.csv", row.names = FALSE)

#########################################################
### by model
library(stringr)
# with 100 trans / MY
dat <- do.call(rbind, df_out)
dat$model <- c("1Char", "2Char", "3Char")[match(str_count(dat$par), unique(str_count(dat$par)))]
# Apply the aggregate function
result <- aggregate(dat[,4:6], 
  by = list(ntips = dat$ntips, 
    type = dat$type, 
    model = dat$model), 
  FUN = function(x) sapply(agg_funcs, function(f) f(x, na.rm = TRUE))
)
result <- cbind(include="all",result)
# Unlist and convert to a data frame for better readability
final_df <- data.frame(result)

# <50 trans / MY
dat_prune_50 <- as.data.frame(dat)[!dat$value > 50,]
# Apply the aggregate function
result <- aggregate(dat_prune_50[,4:6], 
  by = list(ntips = dat_prune_50$ntips, 
    type = dat_prune_50$type, 
    model = dat_prune_50$model), 
  FUN = function(x) sapply(agg_funcs, function(f) f(x, na.rm = TRUE))
)
result <- cbind(include="<50",result)
result$value[7:9,2]/result$value[,2]
final_df <- rbind(final_df, result)

# <10 trans / MY
dat_prune_10 <- as.data.frame(dat)[!dat$value > 10,]
# Apply the aggregate function
result <- aggregate(dat_prune_10[,4:6], 
  by = list(ntips = dat_prune_10$ntips, 
    type = dat_prune_10$type,
    model = dat_prune_10$model), 
  FUN = function(x) sapply(agg_funcs, function(f) f(x, na.rm = TRUE))
)
result <- cbind(include="<10",result)
result$value[7:9,2]/result$value[,2]
final_df <- rbind(final_df, result)

colnames(final_df) <- c("include", "ntips", "type", "model", "estimate", "sim", "diff")

write.csv(final_df, "tables/par-summary-by-type-model.csv", row.names = FALSE)
