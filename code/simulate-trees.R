library(TreeSim)

setwd("~/corhmm-dredge/")

n = c(100, 250, 500)
lambda = 1
mu = 0

for(i in n){
  phy <- sim.bd.taxa(i, 10, lambda, mu, frac = 1, complete = FALSE, stochsampling = FALSE)
  for(j in 1:length(phy)){
    phy[[j]]$edge.length <- phy[[j]]$edge.length/max(branching.times(phy[[j]]))
  }
  file_name <- paste0("trees/tree_", i, ".tre")
  write.tree(phy, file_name)
}
