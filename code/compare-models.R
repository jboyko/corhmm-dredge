# ============================================================================
# 1. CONDITIONAL PERFORMANCE ANALYSIS
# ============================================================================

cat("\n=== CONDITIONAL PERFORMANCE: WHEN ONE MODEL FAILS, HOW DO OTHERS DO? ===\n")

test_cols <- c("test_1", "test_2", "test_3", "test_4")
model_names <- unique(result_df_all$model)

# Function to analyze conditional performance for a specific test
analyze_conditional_performance <- function(data, test_name) {
  
  cat(sprintf("\n--- ANALYSIS FOR %s ---\n", toupper(test_name)))
  
  results <- list()
  
  for(focal_model in model_names) {
    
    cat(sprintf("\nWhen %s gets %s WRONG:\n", focal_model, test_name))
    
    # Get cases where focal model failed this test
    focal_fails <- data$model == focal_model & data[[test_name]] == FALSE
    
    if(sum(focal_fails) == 0) {
      cat(sprintf("  %s never fails %s!\n", focal_model, test_name))
      next
    }
    
    # Find the same test instances (same row indices divided by 5, since we have 5 models per case)
    fail_indices <- which(focal_fails)
    
    # Get corresponding instances for other models
    # Assuming data is structured with models cycling through l0,l1,l2,er,sa repeatedly
    base_cases <- ((fail_indices - 1) %/% 5) * 5  # Get base index for each case
    
    other_model_performance <- data.frame()
    
    for(other_model in model_names) {
      if(other_model == focal_model) next
      
      # Find indices for this other model in the same cases
      other_indices <- base_cases + match(other_model, model_names)
      
      # Remove indices that are out of bounds
      other_indices <- other_indices[other_indices <= nrow(data)]
      
      if(length(other_indices) > 0) {
        success_rate <- mean(data[other_indices, test_name]) * 100
        n_cases <- length(other_indices)
        
        cat(sprintf("  %s succeeds %.1f%% of the time (n=%d)\n", 
          other_model, success_rate, n_cases))
        
        other_model_performance <- rbind(other_model_performance,
          data.frame(
            other_model = other_model,
            success_rate = success_rate,
            n_cases = n_cases
          ))
      }
    }
    
    results[[focal_model]] <- other_model_performance
  }
  
  return(results)
}

# Analyze each test
conditional_results <- list()
for(test in test_cols) {
  conditional_results[[test]] <- analyze_conditional_performance(result_df_all, test)
}

# ============================================================================
# 2. SUMMARY STATISTICS
# ============================================================================

cat("\n=== SUMMARY: MODEL REDUNDANCY ANALYSIS ===\n")

# Create summary tables
for(test in test_cols) {
  cat(sprintf("\n%s REDUNDANCY MATRIX:\n", toupper(test)))
  cat("When model in ROW fails, success rate of model in COLUMN:\n\n")
  
  # Create matrix for this test
  redundancy_matrix <- matrix(NA, length(model_names), length(model_names))
  rownames(redundancy_matrix) <- model_names
  colnames(redundancy_matrix) <- model_names
  
  for(focal_model in model_names) {
    if(!is.null(conditional_results[[test]][[focal_model]])) {
      for(i in 1:nrow(conditional_results[[test]][[focal_model]])) {
        other_model <- conditional_results[[test]][[focal_model]]$other_model[i]
        success_rate <- conditional_results[[test]][[focal_model]]$success_rate[i]
        redundancy_matrix[focal_model, other_model] <- success_rate
      }
    }
  }
  
  # Fill diagonal with NA (can't compare model with itself)
  diag(redundancy_matrix) <- NA
  
  print(round(redundancy_matrix, 1))
}

# ============================================================================
# 3. COMPLEMENTARITY ANALYSIS
# ============================================================================

cat("\n=== MODEL COMPLEMENTARITY ANALYSIS ===\n")

# For each test, find which models are most complementary
for(test in test_cols) {
  cat(sprintf("\n%s - Model Complementarity:\n", toupper(test)))
  
  model_pairs_performance <- data.frame()
  
  for(i in 1:(length(model_names)-1)) {
    for(j in (i+1):length(model_names)) {
      model1 <- model_names[i]
      model2 <- model_names[j]
      
      # Get data for both models
      model1_data <- result_df_all[result_df_all$model == model1, ]
      model2_data <- result_df_all[result_df_all$model == model2, ]
      
      # Ensure we're comparing the same cases
      n_cases <- min(nrow(model1_data), nrow(model2_data))
      if(n_cases > 0) {
        model1_results <- model1_data[1:n_cases, test]
        model2_results <- model2_data[1:n_cases, test]
        
        # Calculate complementarity metrics
        both_succeed <- sum(model1_results & model2_results)
        both_fail <- sum(!model1_results & !model2_results)
        model1_only <- sum(model1_results & !model2_results)
        model2_only <- sum(!model1_results & model2_results)
        
        # At least one succeeds
        at_least_one <- sum(model1_results | model2_results)
        
        model_pairs_performance <- rbind(model_pairs_performance,
          data.frame(
            pair = paste(model1, "+", model2),
            both_succeed = both_succeed,
            both_fail = both_fail,
            model1_only = model1_only,
            model2_only = model2_only,
            at_least_one = at_least_one,
            complementarity = (model1_only + model2_only) / n_cases * 100
          ))
      }
    }
  }
  
  # Sort by complementarity (cases where exactly one succeeds)
  model_pairs_performance <- model_pairs_performance[order(-model_pairs_performance$complementarity), ]
  
  cat("Model pairs ranked by complementarity (% cases where exactly one succeeds):\n")
  print(model_pairs_performance[, c("pair", "complementarity", "at_least_one")])
}

# ============================================================================
# 4. FAILURE PATTERN ANALYSIS
# ============================================================================

cat("\n=== FAILURE PATTERN ANALYSIS ===\n")

for(test in test_cols) {
  cat(sprintf("\n%s - When models fail together vs independently:\n", toupper(test)))
  
  # Reconstruct the data to align models by case
  aligned_data <- data.frame()
  
  # Group by case (assuming 5 models per case)
  for(case_start in seq(1, nrow(result_df_all), by = 5)) {
    case_end <- min(case_start + 4, nrow(result_df_all))
    case_data <- result_df_all[case_start:case_end, ]
    
    if(nrow(case_data) == 5) {  # Complete case
      case_results <- setNames(case_data[[test]], case_data$model)
      aligned_data <- rbind(aligned_data, as.data.frame(t(case_results)))
    }
  }
  
  if(nrow(aligned_data) > 0) {
    # Count failure patterns
    failure_patterns <- data.frame()
    
    for(i in 1:nrow(aligned_data)) {
      failed_models <- names(aligned_data)[!aligned_data[i, ]]
      n_failures <- length(failed_models)
      
      failure_patterns <- rbind(failure_patterns,
        data.frame(
          case = i,
          n_failures = n_failures,
          failed_models = paste(failed_models, collapse = ",")
        ))
    }
    
    # Summarize failure patterns
    failure_summary <- table(failure_patterns$n_failures)
    cat("Number of models failing simultaneously:\n")
    print(failure_summary)
    
    # Most common failure combinations
    if(nrow(failure_patterns) > 0) {
      common_failures <- sort(table(failure_patterns$failed_models), decreasing = TRUE)
      cat("\nMost common failure combinations:\n")
      print(head(common_failures, 10))
    }
  }
}

# ============================================================================
# 5. VISUALIZATIONS
# ============================================================================

# Set up plotting area
par(mfrow = c(2, 2), mar = c(8, 4, 4, 2))

# Plot 1: Complementarity heatmap for test_1
test_name <- "test_1"
if(!is.null(conditional_results[[test_name]])) {
  redundancy_matrix <- matrix(NA, length(model_names), length(model_names))
  rownames(redundancy_matrix) <- model_names
  colnames(redundancy_matrix) <- model_names
  
  for(focal_model in model_names) {
    if(!is.null(conditional_results[[test_name]][[focal_model]])) {
      for(i in 1:nrow(conditional_results[[test_name]][[focal_model]])) {
        other_model <- conditional_results[[test_name]][[focal_model]]$other_model[i]
        success_rate <- conditional_results[[test_name]][[focal_model]]$success_rate[i]
        redundancy_matrix[focal_model, other_model] <- success_rate
      }
    }
  }
  
  # Only plot if we have data
  if(any(!is.na(redundancy_matrix))) {
    image(1:ncol(redundancy_matrix), 1:nrow(redundancy_matrix),
      t(redundancy_matrix[nrow(redundancy_matrix):1,]),
      col = heat.colors(20, rev = TRUE),
      main = paste("When Model Fails", test_name),
      sub = "Success Rate of Other Models (%)",
      xlab = "Other Models", ylab = "Failing Models",
      axes = FALSE)
    axis(1, 1:ncol(redundancy_matrix), colnames(redundancy_matrix), las = 2)
    axis(2, 1:nrow(redundancy_matrix), rev(rownames(redundancy_matrix)), las = 2)
    
    # Add values
    for(i in 1:nrow(redundancy_matrix)) {
      for(j in 1:ncol(redundancy_matrix)) {
        if(!is.na(redundancy_matrix[i,j])) {
          text(j, nrow(redundancy_matrix) + 1 - i, 
            round(redundancy_matrix[i,j], 0), cex = 0.7)
        }
      }
    }
  }
}

# Plot 2-4: Similar heatmaps for other tests
for(plot_num in 2:4) {
  test_name <- test_cols[plot_num]
  
  if(!is.null(conditional_results[[test_name]])) {
    redundancy_matrix <- matrix(NA, length(model_names), length(model_names))
    rownames(redundancy_matrix) <- model_names
    colnames(redundancy_matrix) <- model_names
    
    for(focal_model in model_names) {
      if(!is.null(conditional_results[[test_name]][[focal_model]])) {
        for(i in 1:nrow(conditional_results[[test_name]][[focal_model]])) {
          other_model <- conditional_results[[test_name]][[focal_model]]$other_model[i]
          success_rate <- conditional_results[[test_name]][[focal_model]]$success_rate[i]
          redundancy_matrix[focal_model, other_model] <- success_rate
        }
      }
    }
    
    if(any(!is.na(redundancy_matrix))) {
      image(1:ncol(redundancy_matrix), 1:nrow(redundancy_matrix),
        t(redundancy_matrix[nrow(redundancy_matrix):1,]),
        col = heat.colors(20, rev = TRUE),
        main = paste("When Model Fails", test_name),
        sub = "Success Rate of Other Models (%)",
        xlab = "Other Models", ylab = "Failing Models",
        axes = FALSE)
      axis(1, 1:ncol(redundancy_matrix), colnames(redundancy_matrix), las = 2)
      axis(2, 1:nrow(redundancy_matrix), rev(rownames(redundancy_matrix)), las = 2)
      
      for(i in 1:nrow(redundancy_matrix)) {
        for(j in 1:ncol(redundancy_matrix)) {
          if(!is.na(redundancy_matrix[i,j])) {
            text(j, nrow(redundancy_matrix) + 1 - i, 
              round(redundancy_matrix[i,j], 0), cex = 0.7)
          }
        }
      }
    }
  }
}

# Reset plotting parameters
par(mfrow = c(1, 1), mar = c(5, 4, 4, 2))

# ============================================================================
# 6. KEY INSIGHTS
# ============================================================================

cat("\n=== KEY INSIGHTS ===\n")

cat("This analysis shows:\n")
cat("1. When each model fails a test, how often do the other models succeed?\n")
cat("2. Which model pairs are most complementary (one succeeds when other fails)?\n")
cat("3. Do models tend to fail together or independently?\n")
cat("4. Which models provide the best backup when others fail?\n")

cat("\nHigh values in the heatmaps = good redundancy (other models succeed when one fails)\n")
cat("Low values = models tend to fail on the same cases\n")

cat("\n=== ANALYSIS COMPLETE ===\n")
