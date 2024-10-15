setwd("~/corhmm-dredge/")
raw_dat <- read.csv("tables/raw-parameter-comparison.csv")

head(raw_dat)

raw_dat$ntips <- factor(raw_dat$ntips)
raw_dat$type <- factor(raw_dat$type, levels = c("l0", "l1", "l2", "er"))
levels(raw_dat$type) <- gsub("l0", "Mk", levels(raw_dat$type))


agg_funcs <- list(mean = mean, varriance = var, median = median)
dat <- raw_dat
# dat[,c(4,5)] <- log(dat[,c(4,5)])
# Apply the aggregate function
result <- aggregate(dat[,c(4:5)], 
  by = list(ntips = dat$ntips, type = dat$type), 
  FUN = function(x) sapply(agg_funcs, function(f) f(x, na.rm = TRUE))
)
print(result)

type_colors <- setNames(RColorBrewer::brewer.pal(4, "Set1"), c("Mk", "l1", "l2", "er"))
transp_cols <- sapply(type_colors, adjustcolor, alpha.f = 0.1)

x_index <- c(seq(from = 1, to = 3, length.out = 4), seq(from = 5, to = 7, length.out = 4), seq(from = 9, to = 11, length.out = 4))

pdf("plots/par-plot.pdf", width = 5)
par(mar = c(5, 5, 5, 5))
plot(1, type = "n", xlim = c(0,12),
  ylim = c(-5, 10), xaxt = "n", xlab = "Ntaxa by Regularization Type", ylab = expression(hat(q)[ij]-q[ij]),
  bty = "n")

axis(1, at = x_index[1:4],
  labels = c(NA, NA, NA, NA))
axis(1, at = mean(x_index[1:4]),
  labels = "100", tick = FALSE)

axis(1, at = x_index[5:8],
  labels = c(NA, NA, NA, NA))
axis(1, at = mean(x_index[5:8]),
  labels = "250", tick=FALSE)

axis(1, at = x_index[9:12],
  labels = c(NA, NA, NA, NA))
axis(1, at = mean(x_index[9:12]),
  labels = "500", tick=FALSE)

points(jitter(as.numeric(x_index)[as.numeric(interaction(raw_dat$type, raw_dat$ntips))], amount = 0.1),
  raw_dat$diff,
  pch = 16, col = rep(transp_cols, 3)[as.numeric(interaction(raw_dat$type, raw_dat$ntips))], cex = .5)  # Customize point type, color, and size

boxplot(diff ~ type * ntips, data = raw_dat, 
  at = x_index,
  col = type_colors,
  xlab = "Ntips by Type", ylab = "Diff",
  las = 2, # Rotate axis labels for clarity
  border = "black", 
  ylim=c(0, 10), 
  outline=FALSE,
  add=TRUE, 
  boxwex = 0.5, axes=FALSE)
legend(x=11.5, y=10.2, legend = levels(raw_dat$type), fill = type_colors, border = "black", box.lwd = 1, cex = 0.85, pt.cex = 0.85, xpd = TRUE)
abline(h=0, lty=2)
dev.off()

# plot(1, type = "n", xlim = c(0,12),
#   ylim = c(0, 10), xaxt = "n", xlab = "Ntips by Type", ylab = expression(hat(q)[ij]),
#   bty = "n")
# 
# # points(jitter(c(mean(x_index[1:4]), mean(x_index[5:8]), mean(x_index[9:12]))[as.numeric(dat$ntips)], amount = 1.5), 
# #   dat$true, pch=16, col = adjustcolor("grey", alpha.f = 0.005), cex = 2)
# 
# axis(1, at = x_index[1:4],
#   labels = c(NA, NA, NA, NA))
# axis(1, at = mean(x_index[1:4]),
#   labels = "100", tick = FALSE)
# 
# axis(1, at = x_index[5:8],
#   labels = c(NA, NA, NA, NA))
# axis(1, at = mean(x_index[5:8]),
#   labels = "250", tick=FALSE)
# 
# axis(1, at = x_index[9:12],
#   labels = c(NA, NA, NA, NA))
# axis(1, at = mean(x_index[9:12]),
#   labels = "500", tick=FALSE)
# 
# points(jitter(as.numeric(x_index)[as.numeric(interaction(raw_dat$type, raw_dat$ntips))], amount = 0.1),
#   raw_dat$value,
#   pch = 16, col = rep(transp_cols, 3)[as.numeric(interaction(raw_dat$type, raw_dat$ntips))], cex = .5)  # Customize point type, color, and size
# 
# boxplot(value ~ type * ntips, data = raw_dat, 
#   at = x_index,
#   col = type_colors,
#   xlab = "Ntips by Type", ylab = "Diff",
#   las = 2, # Rotate axis labels for clarity
#   border = "black", 
#   ylim=c(0, 10), 
#   outline=FALSE,
#   add=TRUE, 
#   boxwex = 0.5, axes=FALSE)
# legend("topright", legend = c("unreg", levels(raw_dat$type)[-1]), fill = type_colors, border = "darkblue", bty = "n")




