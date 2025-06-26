break_size <- 1

fit1_data <- dredge_fitsSA$sa_fits[[1]]
scores1 <- fit1_data$every_score[!is.na(fit1_data$every_move)]
moves1 <- fit1_data$every_move[!is.na(fit1_data$every_move)]
accepted1 <- fit1_data$accepted[!is.na(fit1_data$every_move)]
accepted1[is.na(accepted1)] <- 0

fit2_data <- dredge_fitsSA$sa_fits[[2]]
scores2 <- fit2_data$every_score[!is.na(fit2_data$every_move)]
moves2 <- fit2_data$every_move[!is.na(fit2_data$every_move)]
accepted2 <- fit2_data$accepted[!is.na(fit2_data$every_move)]
accepted2[is.na(accepted2)] <- 0
accepted2[1] <- 1

na_gap <- rep(NA, break_size)

combined_scores <- c(scores1, na_gap, scores2)
combined_moves <- c(moves1, na_gap, moves2)
combined_accepted <- c(accepted1, na_gap, accepted2)

len1 <- length(scores1)
total_len <- length(combined_scores)

all_possible_moves <- unique(c(moves1, moves2), na.rm = TRUE)


# move_colors_palette <- c("#2E86AB", "#A23B72", "#F18F01", "#C73E1D", "#6A994E", "#7209B7", "#F77F00", "#003566")
move_colors_palette <- c("grey90", "#2E86AB", "#A23B72", "#7209B7","grey90")
move_colors <- setNames(move_colors_palette[1:length(all_possible_moves)], all_possible_moves)

y_range <- max(combined_scores, na.rm = TRUE) - min(combined_scores, na.rm = TRUE)
y_bottom <- min(combined_scores, na.rm = TRUE) - y_range * 0.25
y_top <- max(combined_scores, na.rm = TRUE) + y_range * 0.05

plot(seq_along(combined_scores), combined_scores, 
  type = "n",
  xlab = "", 
  ylab = "",
  main = "",
  cex.main = 1.4,
  cex.lab = 1.2,
  cex.axis = 1.0,
  font.main = 2,
  las = 1,
  ylim = c(y_bottom, y_top),
  bty = "l",
  xaxs = "i",
  xlim = c(0.5, total_len + 0.5)
)

abline(h = pretty(combined_scores, n = 8), col = "gray90", lty = 1, lwd = 0.5)

abline(v = len1 + break_size / 2, col = "gray40", lty = "dashed", lwd = 1.5)

# mtext("Fit 1", side = 1, line = 2.5, at = len1 / 2, cex = 1.1)
# mtext("Fit 2", side = 1, line = 2.5, at = len1 + break_size + length(scores2)/2, cex = 1.1)

lines(seq_along(combined_scores), combined_scores, lwd = 2, col = "gray50")

accepted_idx <- which(combined_accepted == 1)
rejected_idx <- which(combined_accepted == 0)

points(accepted_idx, combined_scores[accepted_idx], 
  col = "darkgreen", pch = 19, cex = 1.0)

points(rejected_idx, combined_scores[rejected_idx], 
  col = "firebrick", pch = 4, cex = 0.8, lwd = 2)

rect_bottom <- y_bottom + y_range * 0.05
rect_top <- y_bottom + y_range * 0.15

for(i in seq_along(combined_scores)) {
  if(!is.na(combined_moves[i])) {
    rect(xleft = i - 0.5, 
      xright = i + 0.5,
      ybottom = rect_bottom,
      ytop = rect_top,
      col = move_colors[combined_moves[i]],
      border = 'grey20')
  }
}

# legend("topright",
#   inset = c(0.02, 0.02),
#   legend = c("Accepted", "Rejected", names(move_colors)),
#   col = c("darkgreen", "firebrick", move_colors),
#   pch = c(19, 4, rep(15, length(move_colors))),
#   pt.cex = c(1.2, 1.2, rep(2, length(move_colors))),
#   lwd = c(NA, 2, rep(NA, length(move_colors))),
#   bty = "n",
#   cex = 0.9,
#   text.col = "gray20")