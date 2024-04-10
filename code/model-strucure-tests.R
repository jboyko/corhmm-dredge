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
corhmm_fit_l1 <- corHMM:::corHMMDredge(phy = phy, data = cor_data, max.rate.cat = 1, pen_type = "l1", root.p="maddfitz")
print(corhmm_fit_l1$solution)



