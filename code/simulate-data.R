setwd("~/corhmm-dredge/")

library(corHMM)
library(parallel)
library(MASS)
source("code/utils.R")

overwrite = TRUE
nsim = 100

state_table <- read.csv("state-table.csv")
trees <- lapply(dir("trees/", full.names = TRUE), read.tree)

simulation <- state_table[1,1]
nChar <- state_table[1,2]
nStates <- state_table[1,3]
nHidden <- state_table[1,4]

# creates an index mat appropriate for nchar, nstates, and nhidden
index_mat <- get_index_mat(simulation, nChar, nStates, nHidden)

# samples the possible parameter values 
par_table <- get_par_table(index_mat, nSim)

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

res_reg[[74]]
res_unreg[[74]]
rate_mats[[74]]

index_mat$min_df

data_file_name <- paste0("data_", simulation_name, ".RDS")

data_file_found <- data_file_name %in% dir("data/")

# simulate data
if(!data_file_found | overwrite){
  # datasets for each tree
  dat <- lapply(trees, function(x) 
    lapply(1:nsim, function(y) 
      corHMM:::simMarkov(x, Q, root.p)))
  
  # datasets formatted for corhmm
  cor_dat <- lapply(dat, function(x) 
    lapply(x, function(y) 
      data.frame(sp = names(y$TipStates), d = y$TipStates)))
  
  saveRDS(dat, file = paste0("data/", data_file_name))
}

res_file_name <- paste0("res_", simulation_name, ".RDS")

res_file_found <- res_file_name %in% dir("res/")

if(!res_file_found | overwrite){
  # datasets for each tree
  dat <- lapply(trees, function(x) 
    lapply(1:nsim, function(y) 
      corHMM:::simMarkov(x, Q, root.p)))
  
  # datasets formatted for corhmm
  cor_dat <- lapply(dat, function(x) 
    lapply(x, function(y) 
      data.frame(sp = names(y$TipStates), d = y$TipStates)))
  
  saveRDS(dat, file = paste0("data/", data_file_name))
}


res_unreg <- mclapply(cor_dat, function(x) corHMM(phy, x, 1), mc.cores = 10)
res_reg <- mclapply(cor_dat, function(x) corHMM:::corHMMDredge(phy, x, 1, pen_type = "logl1"), mc.cores = 10)

df_unreg <- do.call(rbind, lapply(res_unreg, function(x) c(x$solution)[c(3,2)]))
df_reg <- do.call(rbind, lapply(res_reg, function(x) c(x$solution)[c(3,2)]))

plot_data <- (cbind(df_unreg, df_reg))
colnames(plot_data) <- c("01_unreg", "10_unreg", "01_reg", "10_reg")

bias = colMeans(plot_data) - 1
varr = apply(plot_data, 2, var)
mse = colMeans((plot_data - 1)^2)
rmse = sqrt(colMeans((plot_data - 1)^2))

t(data.frame(bias, varr, mse, rmse))

boxplot(plot_data); abline(h = 1)


library(MCMCpack)
library(nloptr)
library(geiger)
library(TreeSim)
library(parallel)

# Define the prior distributions
q1_prior <- list(rate=10)
q2_prior <- list(rate=10)

# x <- seq(0, 1, length.out=100)
# y <- dexp(x, q1_prior$rate)
# plot(x, y, type="l", xlab="ef", ylab="density")


# Define the posterior distribution function
log_posterior <- function(params, tree, data) {
  q1 <- params[1]
  q2 <- params[2]
  lp_q1 <- dexp(q1, q1_prior$rate, log=TRUE)
  lp_q2 <- dexp(q2, q2_prior$rate, log=TRUE)
  lp_like <- corHMM(tree, data, rate.cat = 1, model = "ARD", p = params, node.states = "none")$loglik
  lp_posterior <- lp_q1 + lp_q2 + lp_like
  if(is.nan(lp_posterior)){
    lp_posterior <- -1e10
  } 
  return(lp_posterior)
}

mcmc_samples <- MCMCmetrop1R(log_posterior, theta.init=c(0.5, 0.5), 
                             mcmc=10000, burnin=1000, thin=10, 
                             verbose=TRUE, tree=phy, data=cor_dat[[84]], logfun=TRUE)
res_unreg_1 <- corHMM(phy, cor_dat[[84]], 1)
res_reg_1 <- corHMM:::corHMMDredge(phy, cor_dat[[84]], 1, pen_type = "logl1")

raftery.diag(mcmc_samples)
plot(mcmc_samples)
summary(mcmc_samples)$quantiles



