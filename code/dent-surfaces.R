library(corHMM)
data(primates)
phy <- multi2di(primates[[1]])
data <- primates[[2]]
corhmm_fit <- corHMM(phy = phy, data = data, rate.cat = 1)

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

corhmm_fit <- corHMM(phy = phy, data = data, rate.cat = 1)
dredge_fit <- corHMM:::corHMMDredge(phy = phy, data = data, 1)

loglik <- corhmm_fit$loglik
best_neglnL <- -loglik
best_par <- corhmm_fit$solution[!is.na(corhmm_fit$solution)]
names(best_par) <- c("rate_01|00", "rate_00|01", "rate_11|01", "rate_01|11")
corhmm_dent <- dent_walk(par=best_par, fn=fn_corHMM, best_neglnL=best_neglnL, nsteps=1000, phy = phy, data = data)
dredge_dent <- dent_walk(par=best_par, fn=fn_corHMM, best_neglnL=best_neglnL, nsteps=1000, phy = phy, data = data)

par(mfrow=c(1,2))
plot(corhmm_dent)
plot(dredge_dent)


