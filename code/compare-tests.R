# ============================================================================
# 1. MODEL PERFORMANCE BY TEST
# ============================================================================

cat("\n=== MODEL PERFORMANCE BY TEST ===\n")

# Calculate success rates for each model on each test
model_performance <- aggregate(. ~ model, 
  data = result_df_all[, c("test_1", "test_2", "test_3", "test_4", "model")], 
  FUN = function(x) mean(x) * 100)

print("Success rates by model and test (%):")
print(model_performance)

# Create a matrix for easier analysis
perf_matrix <- as.matrix(model_performance[, -1])
rownames(perf_matrix) <- model_performance$model
colnames(perf_matrix) <- paste0("test_", 1:4)

print("\nPerformance Matrix (%):")
print(round(perf_matrix, 1))

# ============================================================================
# 2. MODEL RANKINGS BY TEST
# ============================================================================

cat("\n=== MODEL RANKINGS BY TEST ===\n")

# Rank models for each test (1 = best performance)
rankings <- apply(perf_matrix, 2, function(x) rank(-x, ties.method = "average"))

print("Model rankings by test (1 = best):")
print(rankings)

# Average ranking across all tests
avg_ranking <- rowMeans(rankings)
print("\nAverage ranking across all tests:")
print(sort(avg_ranking))

# ============================================================================
# 3. OVERALL MODEL PERFORMANCE
# ============================================================================

cat("\n=== OVERALL MODEL PERFORMANCE ===\n")

# Overall success rate per model (across all tests)
overall_perf <- rowMeans(perf_matrix)
print("Overall performance by model (%):")
print(sort(overall_perf, decreasing = TRUE))

# Count total successes per model
total_tests_per_model <- table(result_df_all$model) * 4  # 4 tests per observation
total_successes <- aggregate(cbind(test_1, test_2, test_3, test_4) ~ model, 
  data = result_df_all, 
  FUN = sum)
total_successes$total <- rowSums(total_successes[, -1])

print("\nTotal successes by model:")
print(total_successes[, c("model", "total")])

# ============================================================================
# 4. MODEL COMPARISON STATISTICS
# ============================================================================

cat("\n=== MODEL COMPARISON STATISTICS ===\n")

# Standard deviation of performance across tests (consistency measure)
model_consistency <- apply(perf_matrix, 1, sd)
print("Performance consistency (lower SD = more consistent):")
print(sort(model_consistency))

# Best and worst test for each model
for(model in rownames(perf_matrix)) {
  best_test <- which.max(perf_matrix[model, ])
  worst_test <- which.min(perf_matrix[model, ])
  
  cat(sprintf("%s: Best on %s (%.1f%%), Worst on %s (%.1f%%)\n",
    model,
    names(best_test), perf_matrix[model, best_test],
    names(worst_test), perf_matrix[model, worst_test]))
}

# ============================================================================
# 5. PAIRWISE MODEL COMPARISONS
# ============================================================================

cat("\n=== PAIRWISE MODEL COMPARISONS ===\n")

# Compare each pair of models across all tests
model_names <- unique(result_df_all$model)
comparison_results <- data.frame()

for(i in 1:(length(model_names)-1)) {
  for(j in (i+1):length(model_names)) {
    model1 <- model_names[i]
    model2 <- model_names[j]
    
    # Count wins for each model across tests
    wins_model1 <- sum(perf_matrix[model1, ] > perf_matrix[model2, ])
    wins_model2 <- sum(perf_matrix[model2, ] > perf_matrix[model1, ])
    ties <- sum(perf_matrix[model1, ] == perf_matrix[model2, ])
    
    comparison_results <- rbind(comparison_results, 
      data.frame(
        comparison = paste(model1, "vs", model2),
        model1_wins = wins_model1,
        model2_wins = wins_model2,
        ties = ties,
        model1_avg = mean(perf_matrix[model1, ]),
        model2_avg = mean(perf_matrix[model2, ])
      ))
  }
}

print("Head-to-head comparisons (wins across 4 tests):")
print(comparison_results)

# ============================================================================
# 6. VISUALIZATIONS
# ============================================================================

# Set up plotting area
par(mfrow = c(2, 2), mar = c(5, 4, 4, 2))

# Plot 1: Performance by Model and Test
barplot(t(perf_matrix), 
  beside = TRUE,
  main = "Model Performance by Test",
  ylab = "Success Rate (%)",
  col = rainbow(4),
  legend.text = colnames(perf_matrix),
  args.legend = list(x = "topright", cex = 0.8),
  ylim = c(0, 100))

# Plot 2: Overall Model Performance
barplot(sort(overall_perf, decreasing = TRUE),
  main = "Overall Model Performance",
  ylab = "Average Success Rate (%)",
  col = "lightblue",
  ylim = c(0, 100))
# Add percentage labels
text(x = 1:length(overall_perf) * 1.2 - 0.5,
  y = sort(overall_perf, decreasing = TRUE) + 2,
  labels = paste0(round(sort(overall_perf, decreasing = TRUE), 1), "%"),
  cex = 0.8)

# Plot 3: Model Consistency (SD)
barplot(sort(model_consistency),
  main = "Model Consistency Across Tests",
  ylab = "Standard Deviation (%)",
  col = "lightcoral",
  sub = "Lower = More Consistent")

# Plot 4: Heatmap of Model Performance
image(1:ncol(perf_matrix), 1:nrow(perf_matrix),
  t(perf_matrix),
  col = heat.colors(20),
  main = "Model Performance Heatmap",
  xlab = "Tests", ylab = "Models",
  axes = FALSE)
axis(1, 1:ncol(perf_matrix), colnames(perf_matrix))
axis(2, 1:nrow(perf_matrix), rownames(perf_matrix))

# Add performance values to heatmap
for(i in 1:nrow(perf_matrix)) {
  for(j in 1:ncol(perf_matrix)) {
    text(j, i, round(perf_matrix[i,j], 0), cex = 0.8)
  }
}

# Reset plotting parameters
par(mfrow = c(1, 1), mar = c(5, 4, 4, 2))

# ============================================================================
# 7. DETAILED MODEL PROFILES
# ============================================================================

cat("\n=== DETAILED MODEL PROFILES ===\n")

for(model in model_names) {
  cat(sprintf("\n--- %s MODEL PROFILE ---\n", toupper(model)))
  
  # Performance on each test
  perf <- perf_matrix[model, ]
  cat("Test Performance:\n")
  for(test in names(perf)) {
    cat(sprintf("  %s: %.1f%%\n", test, perf[test]))
  }
  
  # Overall stats
  cat(sprintf("Overall Average: %.1f%%\n", mean(perf)))
  cat(sprintf("Best Test: %s (%.1f%%)\n", names(which.max(perf)), max(perf)))
  cat(sprintf("Worst Test: %s (%.1f%%)\n", names(which.min(perf)), min(perf)))
  cat(sprintf("Range: %.1f percentage points\n", max(perf) - min(perf)))
  cat(sprintf("Consistency (SD): %.1f\n", sd(perf)))
  cat(sprintf("Average Ranking: %.1f\n", avg_ranking[model]))
}

# ============================================================================
# 8. SUMMARY INSIGHTS
# ============================================================================

cat("\n=== KEY INSIGHTS ===\n")

best_overall <- names(which.max(overall_perf))
worst_overall <- names(which.min(overall_perf))

cat(sprintf("Best overall model: %s (%.1f%% average)\n", 
  best_overall, overall_perf[best_overall]))
cat(sprintf("Worst overall model: %s (%.1f%% average)\n", 
  worst_overall, overall_perf[worst_overall]))

most_consistent <- names(which.min(model_consistency))
least_consistent <- names(which.max(model_consistency))

cat(sprintf("Most consistent model: %s (SD = %.1f)\n", 
  most_consistent, model_consistency[most_consistent]))
cat(sprintf("Least consistent model: %s (SD = %.1f)\n", 
  least_consistent, model_consistency[least_consistent]))

# Find which test is hardest/easiest overall
test_difficulty <- colMeans(perf_matrix)
hardest_test <- names(which.min(test_difficulty))
easiest_test <- names(which.max(test_difficulty))

cat(sprintf("Hardest test overall: %s (%.1f%% average)\n", 
  hardest_test, test_difficulty[hardest_test]))
cat(sprintf("Easiest test overall: %s (%.1f%% average)\n", 
  easiest_test, test_difficulty[easiest_test]))

# Model dominance
cat("\nModel dominance (how often each model is best on each test):\n")
dominance <- apply(perf_matrix, 2, function(x) names(which.max(x)))
for(test in names(dominance)) {
  cat(sprintf("  %s: %s is best\n", test, dominance[test]))
}

cat("\n=== ANALYSIS COMPLETE ===\n")

