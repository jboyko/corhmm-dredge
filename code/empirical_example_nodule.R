library(corHMM)
setwd("~/corhmm-dredge/")

dat <- read.csv("data/empirical_examples/nodules/41467_2024_48036_MOESM4_ESM.csv")
cor_dat <- data.frame(sp = dat$Taxon.in.Tree, nods = dat$State)

phy <- read.tree("data/empirical_examples/nodules/41467_2024_48036_MOESM6_ESM.txt")

rownames(cor_dat) <- cor_dat$sp
cor_dat$sp[(duplicated(cor_dat$sp))]
length(unique(phy$tip.label))

strsplit(phy$tip.label, "\\.")
  