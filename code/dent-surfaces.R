library(corHMM)
data(primates)
phy <- multi2di(primates[[1]])
data <- primates[[2]]

library(dentist)

fn_corHMM <- function(par, phy, data, rate.cat){
  corhmm_fit <- corHMM(phy = phy, data = data, rate.cat = 1, p = par)
  loglik <- corhmm_fit$loglik
  neg_loglik <- -loglik
  return(neg_loglik)
}

fn_Dredge <- function(par, phy, data, rate.cat){
  corhmm_fit <- corHMM(phy = phy, data = data, rate.cat = 1, p = par)
  pen <- sum(log(par))
  loglik <- corhmm_fit$loglik + pen
  neg_loglik <- -loglik
  return(neg_loglik)
}

phy$edge.length <- phy$edge.length + 2e-5

corhmm_fit <- corHMM(phy = phy, data = data, rate.cat = 1, root.p="maddfitz")
dredge_fit_4 <- corHMM:::corHMMDredge(phy = phy, data = data, 4, nstarts = 10, n.cores = 10, lower.bound = .Machine$double.eps)
dredge_fit_3 <- corHMM:::corHMMDredge(phy = phy, data = data, 3, nstarts = 10, n.cores = 10, lower.bound = .Machine$double.eps)
dredge_fit_2 <- corHMM:::corHMMDredge(phy = phy, data = data, 2, nstarts = 10, n.cores = 10, lower.bound = .Machine$double.eps)
dredge_fit_1 <- corHMM:::corHMMDredge(phy = phy, data = data, 1, nstarts = 10, n.cores = 10, lower.bound = .Machine$double.eps)

par(mar=c(.1,.1,.1,.1), mfrow=c(1,3))
plotRECON(dredge_fit_2$phy, dredge_fit_1$states, pie.cex = 1, show.tip.label = FALSE)
plotRECON(dredge_fit_2$phy, dredge_fit_2$states, pie.cex = 1, show.tip.label = FALSE)
plotRECON(dredge_fit_2$phy, dredge_fit_3$states, pie.cex = 1, show.tip.label = FALSE)
plotRECON(dredge_fit_2$phy, dredge_fit_4$states, pie.cex = 1, show.tip.label = FALSE)

cbind(
  d1 = apply(dredge_fit_1$states, 1, which.max),
  d2 = apply(dredge_fit_2$states, 1, which.max),
  d3 = apply(dredge_fit_3$states, 1, which.max),
  d4 = apply(dredge_fit_3$states, 1, which.max)
)

corHMM:::getModelTable(list(dredge_fit_1, dredge_fit_2, dredge_fit_3, dredge_fit_4))

apply(dredge_fit_2$states, 1, which.max)
apply(dredge_fit_3$states, 1, which.max)


dredge_1_profile <- corHMM:::get_batch_profile_lik(corhmm_obj = dredge_fit_1,
                                                 range_factor = 1000,
                                                 n_points = 20,
                                                 ncores = 10,
                                                 dredge = TRUE)
corHMM:::plot_batch_profile_lik(dredge_1_profile, ylim = c(-100, -41))

dredge_2_profile <- corHMM:::get_batch_profile_lik(corhmm_obj = dredge_fit_2,
                                                   range_factor = 1000,
                                                   n_points = 20,
                                                   ncores = 10,
                                                   dredge = TRUE)
corHMM:::plot_batch_profile_lik(dredge_2_profile, ylim = c(-50, -38))

dredge_3_profile <- corHMM:::get_batch_profile_lik(corhmm_obj = dredge_fit_3,
                                                   range_factor = 1000,
                                                   n_points = 20,
                                                   ncores = 10,
                                                   dredge = TRUE)
corHMM:::plot_batch_profile_lik(dredge_3_profile, ylim = c(-60, -37))


trees <- makeSimmap(tree = dredge_fit_2$phy, dredge_fit_2$data, dredge_fit_1$solution, 1, "maddfitz", nSim = 100)
phytools:::describe.simmap(trees)

debug(corHMM:::corHMMDredge)

loglik <- corhmm_fit$loglik
best_neglnL <- -loglik
best_par <- corhmm_fit$solution[!is.na(corhmm_fit$solution)]
names(best_par) <- c("rate_01|00", "rate_00|01", "rate_11|01", "rate_01|11")
corhmm_dent <- dent_walk(par=best_par, fn=fn_corHMM, best_neglnL=best_neglnL, nsteps=1000, phy = phy, data = data)
dredge_dent <- dent_walk(par=best_par, fn=fn_corHMM, best_neglnL=best_neglnL, nsteps=1000, phy = phy, data = data)

par(mfrow=c(1,2))
plot(corhmm_dent)
plot(dredge_dent)


