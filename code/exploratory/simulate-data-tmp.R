library(corHMM)
library(parallel)

nsim = 100
Q <- matrix(c(-1,1,1,-1), 2, 2, byrow = TRUE)
root.p <- c(1,0)

trees <- lapply(dir("trees/", full.names = TRUE), read.tree)

dat <- lapply(1:nsim, function(x) corHMM:::simMarkov(phy, Q, root.p))
cor_dat <- lapply(dat, function(x) data.frame(sp = names(x$TipStates), d = x$TipStates))

# plot(phy, show.tip.label = FALSE)
# tiplabels(pch = 16, col = dat$TipStates, cex = 0.5)

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



