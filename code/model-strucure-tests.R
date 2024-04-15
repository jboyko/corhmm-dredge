library(corHMM)
library(geiger)
library(MASS)

setwd("corhmm-dredge/")
source("code/utils.R")

phy <- sim.bdtree()

index_mat <- get_index_mat(nChar=1, nStates=2, nRateClass=1)
rate_mat <- get_rate_mat(index_mat, c(.1, .5))
sim_data <- get_sim_data(phy, rate_mat, index_mat)
cor_data <- get_formatted_data(sim_data, index_mat)

# LRT
corhmm_fit_ard <- corHMM(phy = phy, data = cor_data, rate.cat = 1, root.p="maddfitz")
corhmm_fit_er <- corHMM(phy = phy, data = cor_data, rate.cat = 1, root.p="maddfitz", model="ER")

LR = -2 * (corhmm_fit_er$loglik-corhmm_fit_ard$loglik)
DF = max(corhmm_fit_ard$index.mat, na.rm=TRUE) - max(corhmm_fit_er$index.mat, na.rm=TRUE)
p_value = pchisq(LR, DF, lower.tail = FALSE)
print(p_value)

# model averaging
mod_table <- corHMM:::getModelTable(list(corhmm_fit_er, corhmm_fit_ard))
avg_solution <- (corhmm_fit_er$solution * mod_table$AICwt[1]) + (corhmm_fit_ard$solution * mod_table$AICwt[2])
print(avg_solution)

# dredging
corhmm_fit_l1 <- corHMM:::corHMMDredge(phy = phy, data = cor_data, max.rate.cat = 1, pen_type = "l1", root.p="maddfitz", lambda = 1)
print(corhmm_fit_l1$solution)

LOO_validation <- function(corhmm_fit, lambda){
  dat <- corhmm_fit$data
  to_rm <- sample(dim(dat)[1], round(dim(dat)[1] * 0.1))
  fit_cor_loo <- function(to_rm_i){
    dat[to_rm_i, 2] <- "?"
    fit <- corHMM:::corHMMDredge(phy = phy, data = dat, max.rate.cat = 1, pen_type = "l1", root.p="maddfitz", lambda = lambda, get.tip.states = TRUE)
    return(fit)
  }
  res <- lapply(to_rm, fit_cor_loo)
  names(res) <- to_rm
  return(res)
}

test <- LOO_validation(corhmm_fit_l1, 1)
fin_table <- cor_data[as.numeric(names(test)), ]
tmp <- t(sapply(1:10, function(x) test[[x]]$tip.states[as.numeric(names(test))[x], ]))
L1 <- cbind(fin_table, tmp)
test <- LOO_validation(corhmm_fit_l1, 0.5)
fin_table <- cor_data[as.numeric(names(test)), ]
tmp <- t(sapply(1:10, function(x) test[[x]]$tip.states[as.numeric(names(test))[x], ]))
L0.5 <- cbind(fin_table, tmp)
test <- LOO_validation(corhmm_fit_l1, 0)
fin_table <- cor_data[as.numeric(names(test)), ]
tmp <- t(sapply(1:10, function(x) test[[x]]$tip.states[as.numeric(names(test))[x], ]))
L0 <- cbind(fin_table, tmp)



