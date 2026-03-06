# ==== Load Libraries ====
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggrepel)


output_dir <- "./inferencePlot/"

# ==== Create Output Folder if Not Exists ====
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# ==== Analysis Function ====
analyse_tree <- function(model_name, tree, param) {
  required_param_cols <- c("r_value", "minD_value", "minVar_value")
  if (!all(required_param_cols %in% names(param))) {
    stop("Missing param columns: expected r_value, minD_value, minVar_value")
  }
  
  param_combination <- param[, required_param_cols, drop = FALSE]
  names(param_combination) <- c("r_Value", "minD", "minVar")
  param_combination <- param_combination %>%
    mutate(across(everything(), as.numeric))
  if (!"r_Value" %in% names(tree)) {
    stop("Missing r_Value column in tree data")
  }
  tree$r_Value <- as.numeric(as.character(tree$r_Value))
  
  tree_has_minD <- "minD" %in% names(tree) && any(!is.na(tree$minD) & tree$minD != "")
  tree_has_minVar <- "minVar" %in% names(tree) && any(!is.na(tree$minVar) & tree$minVar != "")
  
  if (!tree_has_minD || !tree_has_minVar) {
    merged_tree <- merge(tree, param_combination, by = "r_Value", all.x = TRUE)
    merged_tree <- merged_tree %>%
      rename(minD = minD.y, minVar = minVar.y) %>%
      select(-minD.x, -minVar.x, everything())
  } else {
    tree <- tree %>%
      mutate(across(c(minD, minVar), as.numeric))
    merged_tree <- merge(tree, param_combination, by = c("r_Value", "minD", "minVar"), all.x = TRUE)
  }
  
  # Add model name based on script
  merged_tree$model <- case_when(
    grepl("TruewtErrModel", merged_tree$script) ~ "TruewtErrModel",
    grepl("TrueErrModel", merged_tree$script) ~ "TrueErrModel",
    grepl("wtErrModel", merged_tree$script) ~ "wtErrModel",
    grepl("ErrModel", merged_tree$script) ~ "ErrModel",
    TRUE ~ "Other"
  )
  
  merged_tree$param_combo <- paste0("D", merged_tree$minD, "_V", merged_tree$minVar)
  
  # Summary data by param combination
  summary_data <- merged_tree %>%
    group_by(model, param_combo) %>%
    summarise(
      mean_RF_dist = mean(RF_Dist, na.rm = TRUE),
      mean_height_error = mean(MAE_treeH, na.rm = TRUE),
      mean_length_error = mean(MAE_treeL, na.rm = TRUE),
      .groups = 'drop')
  
  # Mean RF distance by parameter combinations (points + lines)
  true_avg <- merged_tree %>%
  filter(grepl("^True", model)) %>%      
  group_by(model) %>%
  summarise(mean_RF = mean(RF_Dist, na.rm = TRUE), .groups = "drop")
    x_right <- length(unique(merged_tree$param_combo)) + 1

    p4_mean <- ggplot() +
      geom_point(
        data = merged_tree,
        aes(x = param_combo, y = RF_Dist, color = model),
        position = position_jitter(width = 0.2, height = 0, seed = 1),
        alpha = 0.4, size = 1.5
      ) +
      geom_point(
        data = summary_data,
        aes(x = param_combo, y = mean_RF_dist, color = model, group = model),
        size = 3
      ) +
      geom_hline(
        data = true_avg,
        aes(yintercept = mean_RF, color = model),
        linewidth = 0.9,
        alpha = 0.9
      ) +
      geom_text_repel(
          data = true_avg,
          aes(
            x = x_right,
            y = mean_RF,
            label = paste0(model, ": ", round(mean_RF, 3)),
            color = model
          ),
          direction = "y",
          hjust = 0,
          size = 4,
          segment.color = NA,
          inherit.aes = FALSE
        )+
      geom_point(
        data = true_avg,
        aes(x = Inf, y = mean_RF, color = model),
        size = 3,
        inherit.aes = FALSE
      ) +
      coord_cartesian(clip = "off") +
      theme_minimal() +
      labs(
        title = paste(model_name, "- RF Distance by Parameter Combination"),
        x = "Parameter Combination (minD_minVar)",
        y = "RF Distance"
      ) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.margin = margin(5.5, 25, 5.5, 5.5) 
      ) + scale_x_discrete(expand = expansion(mult = c(0.05, 0.3))) + 
      scale_color_brewer(type = "qual", palette = "Set1")
  
  ggsave(filename = paste0(output_dir, model_name, "_RF_by_Params_withRaw.png"), 
         plot = p4_mean, width = 12, height = 8, dpi = 300)
  
  
  # Mean height relative error by parameter combinations
  p5_mean <- ggplot() +
    geom_point(data = merged_tree, 
               aes(x = param_combo, y = MAE_treeH, color = model), 
               position = position_jitter(width = 0.2, height = 0, seed = 1),
               alpha = 0.4, size = 1.5) +
    geom_point(data = summary_data, 
               aes(x = param_combo, y = mean_height_error, color = model, group = model),
               size = 3) +
    geom_line(data = summary_data, 
              aes(x = param_combo, y = mean_height_error, color = model, group = model)) +
    theme_minimal() +
    labs(title = paste(model_name, "- Tree Height Error by Parameter Combination"), 
         x = "Parameter Combination", 
         y = "Height Error") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_color_brewer(type = "qual", palette = "Set1")
  
  ggsave(filename = paste0(output_dir, model_name, "_Height_Error_by_Params_withRaw.png"), 
         plot = p5_mean, width = 12, height = 8, dpi = 300)
  
  
  # Mean length relative error by parameter combinations
  p6_mean <- ggplot() +
    geom_point(data = merged_tree, 
               aes(x = param_combo, y = MAE_treeL, color = model), 
               position = position_jitter(width = 0.2, height = 0, seed = 1),
               alpha = 0.4, size = 1.5) +
    geom_point(data = summary_data, 
               aes(x = param_combo, y = mean_length_error, color = model, group = model),
               size = 3) +
    geom_line(data = summary_data, 
              aes(x = param_combo, y = mean_length_error, color = model, group = model)) +
    theme_minimal() +
    labs(title = paste(model_name, "- Tree Length Error by Parameter Combination"), 
         x = "Parameter Combination", 
         y = "Length Error") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_color_brewer(type = "qual", palette = "Set1")
  
  ggsave(filename = paste0(output_dir, model_name, "_Length_Error_by_Params_withRaw.png"), 
         plot = p6_mean, width = 12, height = 8, dpi = 300)
  
  
  # ============================
  # === Combined Accuracy Plot (boxplot, keep original) ===
  # ============================
  
  accuracy_long <- merged_tree %>%
    select(model, param_combo, RF_Dist, MAE_treeH, MAE_treeL) %>%
    pivot_longer(cols = c(RF_Dist, MAE_treeH, MAE_treeL), 
                 names_to = "Metric", values_to = "Value")
  
  p7 <- ggplot(accuracy_long, aes(x = param_combo, y = Value, fill = model)) +
    geom_boxplot() +
    theme_minimal() +
    labs(title = paste(model_name, "- Combined Accuracy Metrics"), 
         x = "Parameter Combination", y = "Metric Value") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    facet_grid(Metric ~ ., scales = "free_y") +
    scale_fill_brewer(type = "qual", palette = "Set1")
  
  ggsave(filename = paste0(output_dir, model_name, "_Combined_Accuracy.png"), 
         plot = p7, width = 14, height = 10, dpi = 300)
  
  # ============================
  # === Summary Plots ===
  # ============================
  
  # Plot mean RF distance by model and parameters (line + point)
  # Make sure required variables are numeric
  merged_tree <- merged_tree %>%
    mutate(across(c(minD, minVar, RF_Dist), as.numeric))
  
  # Subset data (optional if needed to filter models or remove NA)
  rf_data <- merged_tree %>%
    filter(model %in% c("ErrModel", "wtErrModel"),
           !is.na(RF_Dist), !is.na(minD), !is.na(minVar))
  
  # Fit the regression model
  lm_rf <- lm(RF_Dist ~ minD + minVar, data = rf_data)
  
  # Predict values from the model
  rf_data$predicted <- predict(lm_rf, newdata = rf_data)
  
  # Order param_combo based on minD then minVar
  rf_data$param_combo_ordered <- factor(
    rf_data$param_combo,
    levels = rf_data %>% arrange(minD, minVar) %>% pull(param_combo) %>% unique()
  )
  
  # Create regression summary text
  reg_summary <- paste0(
    "y = ", round(coef(lm_rf)[1], 3), " + ",
    round(coef(lm_rf)[2], 3), "·minD + ",
    round(coef(lm_rf)[3], 3), "·minVar\n",
    "R² = ", round(summary(lm_rf)$r.squared, 3)
  )
  
  # Plot
  p_rf <- ggplot(rf_data, aes(x = param_combo_ordered)) +
    geom_point(aes(y = RF_Dist, color = model), size = 2, 
               position = position_jitter(width = 0.3), alpha = 0.7) +
    geom_line(aes(y = predicted, group = 1), color = "black", size = 1) +
    annotate("text", x = Inf, y = Inf, label = reg_summary,
             hjust = 1.1, vjust = 1.1, size = 4, color = "black") +
    theme_minimal() +
    labs(
      title = "Normalised RF Distance (Scale 0 to 1)",
      subtitle = paste0("Ordered by minD then minVar for ", model_name, " model. Regression: RF_Dist ~ minD + minVar"),
      x = "Parameter Combination (minD_minVar)",
      y = "Normalised RF Distance"
    ) +
    scale_color_brewer(palette = "Set2") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Save
  ggsave(paste0(output_dir, model_name, "_Mean_RF_Performance.png"),
         plot = p_rf, width = 12, height = 6, dpi = 300)

  
  # # Plot mean accuracy errors by model and parameters
  # error_long <- summary_data %>%
  #   pivot_longer(cols = c(mean_height_error, mean_length_error), 
  #                names_to = "Error_Type", values_to = "Mean_Error")
  # 
  # p9 <- ggplot(error_long, aes(x = param_combo, y = Mean_Error, color = model, group = model)) +
  #   geom_point(size = 3) +
  #   geom_line() +
  #   theme_minimal() +
  #   labs(title = paste(model_name, "- Mean Estimation Errors by Model and Parameters"), 
  #        x = "Parameter Combination", y = "Mean Relative Error (%)",
  #        subtitle = "Lower values indicate better accuracy") +
  #   theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  #   facet_wrap(~Error_Type, scales = "free_y", 
  #              labeller = labeller(Error_Type = c("mean_height_error" = "Height Error", 
  #                                                 "mean_length_error" = "Length Error"))) +
  #   scale_color_brewer(type = "qual", palette = "Set1")
  # 
  # ggsave(filename = paste0(output_dir, model_name, "_Mean_Accuracy_Performance.png"), 
  #        plot = p9, width = 14, height = 8, dpi = 300)
  # 
  # cat("Generated", 11, "plots for", model_name, "model in", output_dir, "\n")
  # 
  
  # # === Misspecification Analysis ===
  # err_data <- merged_tree %>%
  #   filter(model %in% c("TrueErrModel", "ErrModel")) %>%
  #   select(r_Value, param_combo, model, minD, minVar, RF_Dist, MAE_treeH, MAE_treeL)
  # 
  # wt_data <- merged_tree %>%
  #   filter(model %in% c("TruewtErrModel", "wtErrModel")) %>%
  #   select(r_Value, param_combo, model, minD, minVar, RF_Dist, MAE_treeH, MAE_treeL)
  # 
  # err_wide <- pivot_wider(err_data, id_cols = c(r_Value, param_combo, minD, minVar),
  #                         names_from = model, values_from = c(RF_Dist, MAE_treeH, MAE_treeL))
  # 
  # wt_wide <- pivot_wider(wt_data, id_cols = c(r_Value, param_combo, minD, minVar),
  #                        names_from = model, values_from = c(RF_Dist, MAE_treeH, MAE_treeL))
  # 
  # misspec_r_level <- bind_rows(
  #   err_wide %>% mutate(
  #     Model = "ErrModel",
  #     RF_Diff = abs(RF_Dist_ErrModel - RF_Dist_TrueErrModel),
  #     MAE_treeH = abs(MAE_treeH_ErrModel - MAE_treeH_TrueErrModel),
  #     MAE_treeL = abs(MAE_treeL_ErrModel - MAE_treeL_TrueErrModel)
  #   ) %>% select(r_Value, param_combo, minD, minVar, Model, RF_Diff, MAE_treeH, MAE_treeL),
  #   wt_wide %>% mutate(
  #     Model = "wtErrModel",
  #     RF_Diff = abs(RF_Dist_wtErrModel - RF_Dist_TruewtErrModel),
  #     MAE_treeH = abs(MAE_treeH_wtErrModel - MAE_treeH_TruewtErrModel),
  #     MAE_treeL = abs(MAE_treeL_wtErrModel - MAE_treeL_TruewtErrModel)
  #   ) %>% select(r_Value, param_combo, minD, minVar, Model, RF_Diff, MAE_treeH, MAE_treeL)
  # )
  # 
  # misspec_r_long <- pivot_longer(misspec_r_level,
  #                                cols = c(RF_Diff, MAE_treeH, MAE_treeL),
  #                                names_to = "Error_Type",
  #                                values_to = "Relative") %>%
  #   mutate(Error_Type = recode(Error_Type,
  #                              "RF_Diff" = "RF Distance",
  #                              "MAE_treeH" = "Height Error",
  #                              "MAE_treeL" = "Length Error"))
  # 
  # # ==== RF Distance Regression ====
  # rf_data <- misspec_r_long %>%
  #   filter(Error_Type == "RF Distance", !is.na(AbsDifference), !is.na(minD), !is.na(minVar))
  # 
  # # Fit one linear model for RF distance misspecification
  # lm_rf <- lm(AbsDifference ~ minD + minVar, data = rf_data)
  # rf_data$predicted <- predict(lm_rf, newdata = rf_data)
  # 
  # # Order param_combo by minD then minVar
  # rf_data$param_combo_ordered <- factor(
  #   rf_data$param_combo,
  #   levels = rf_data %>% arrange(minD, minVar) %>% pull(param_combo) %>% unique()
  # )
  # 
  # # Regression summary annotation
  # reg_summary <- paste0(
  #   "y = ", round(coef(lm_rf)[1], 3), " + ",
  #   round(coef(lm_rf)[2], 3), "·minD + ",
  #   round(coef(lm_rf)[3], 3), "·minVar\n",
  #   "R² = ", round(summary(lm_rf)$r.squared, 3)
  # )
  # 
  # # Plot
  # p_rf <- ggplot(rf_data, aes(x = param_combo_ordered)) +
  #   #geom_smooth(aes(y = AbsDifference, group = 1), method = "loess", se = FALSE, color = "blue")+ 
  #   geom_point(aes(y = AbsDifference, color = Model), size = 2, position = position_jitter(width = 0.3), alpha = 0.7) +
  #   geom_line(aes(y = predicted, group = 1), color = "black", size = 1) +
  #   annotate("text", x = Inf, y = Inf, label = reg_summary,
  #            hjust = 1.1, vjust = 1.1, size = 4, color = "black") +
  #   theme_minimal() +
  #   labs(
  #     title = "Relative RF Distance (called alignment - true alignment)",
  #     subtitle = paste0("Ordered by minD then minVar for ", model_name, " model. Regression: Relative RF Distance ~ minD + minVar"),
  #     x = "Parameter Combination",
  #     y = "Relative RF Distance"
  #   ) +
  #   scale_color_brewer(palette = "Set2") +
  #   theme(axis.text.x = element_text(angle = 45, hjust = 1))
  # 
  # ggsave(paste0(output_dir, model_name, "_relative_RF_Distance_withRegression.png"),
  #        p_rf, width = 12, height = 6, dpi = 300)
  # 
  # 
  # # === Regression Plots (Height & Length Error) ===
  # # ==== Height Error Regression ====
  # height_data <- misspec_r_long %>%
  #   filter(Error_Type == "Height Error", !is.na(AbsDifference), !is.na(minD), !is.na(minVar))
  # 
  # # Fit global linear model for all models combined
  # lm_height <- lm(AbsDifference ~ minD + minVar, data = height_data)
  # height_data$predicted <- predict(lm_height, newdata = height_data)
  # 
  # # Order x-axis
  # height_data$param_combo_ordered <- factor(
  #   height_data$param_combo,
  #   levels = height_data %>% arrange(minD, minVar) %>% pull(param_combo) %>% unique()
  # )
  # 
  # # Regression annotation
  # reg_summary_height <- paste0(
  #   "y = ", round(coef(lm_height)[1], 3), " + ",
  #   round(coef(lm_height)[2], 3), "·minD + ",
  #   round(coef(lm_height)[3], 3), "·minVar\n",
  #   "R² = ", round(summary(lm_height)$r.squared, 3)
  # )
  # 
  # Plot
  # p_height_reg <- ggplot(height_data, aes(x = param_combo_ordered)) +
  #   geom_point(aes(y = AbsDifference, color = Model), size = 2, position = position_jitter(width = 0.3), alpha = 0.7) +
  #   geom_line(aes(y = predicted, group = 1), color = "black", size = 1) +
  #   annotate("text", x = Inf, y = Inf, label = reg_summary_height,
  #            hjust = 1.1, vjust = 1.1, size = 4, color = "black") +
  #   theme_minimal() +
  #   labs(
  #     title = paste(model_name, "- Tree Height Misspecification"),
  #     subtitle = "Ordered by minD then minVar | Regression: AbsDifference ~ minD + minVar",
  #     x = "Parameter Combination",
  #     y = "Absolute Height Error Difference (%)"
  #   ) +
  #   scale_color_brewer(palette = "Set2") +
  #   theme(axis.text.x = element_text(angle = 45, hjust = 1))
  # 
  # ggsave(paste0(output_dir, model_name, "_Misspec_Height_Error_perR_withRegression.png"),
  #        p_height_reg, width = 12, height = 6, dpi = 300)
  # 
  # 
  # # ==== Length Error Regression ====
  # length_data <- misspec_r_long %>%
  #   filter(Error_Type == "Length Error", !is.na(AbsDifference), !is.na(minD), !is.na(minVar))
  # 
  # # Fit global linear model
  # lm_length <- lm(AbsDifference ~ minD + minVar, data = length_data)
  # length_data$predicted <- predict(lm_length, newdata = length_data)
  # 
  # # Order x-axis
  # length_data$param_combo_ordered <- factor(
  #   length_data$param_combo,
  #   levels = length_data %>% arrange(minD, minVar) %>% pull(param_combo) %>% unique()
  # )
  # 
  # # Regression annotation
  # reg_summary_length <- paste0(
  #   "y = ", round(coef(lm_length)[1], 3), " + ",
  #   round(coef(lm_length)[2], 3), "·minD + ",
  #   round(coef(lm_length)[3], 3), "·minVar\n",
  #   "R² = ", round(summary(lm_length)$r.squared, 3)
  # )
  # 
  # # Plot
  # p_length_reg <- ggplot(length_data, aes(x = param_combo_ordered)) +
  #   geom_point(aes(y = AbsDifference, color = Model), size = 2, position = position_jitter(width = 0.3), alpha = 0.7) +
  #   geom_line(aes(y = predicted, group = 1), color = "black", size = 1) +
  #   annotate("text", x = Inf, y = Inf, label = reg_summary_length,
  #            hjust = 1.1, vjust = 1.1, size = 4, color = "black") +
  #   theme_minimal() +
  #   labs(
  #     title = paste(model_name, "- Tree Length Misspecification"),
  #     subtitle = "Ordered by minD then minVar | Regression: AbsDifference ~ minD + minVar",
  #     x = "Parameter Combination",
  #     y = "Absolute Length Error Difference (%)"
  #   ) +
  #   scale_color_brewer(palette = "Set2") +
  #   theme(axis.text.x = element_text(angle = 45, hjust = 1))
  # 
  # ggsave(paste0(output_dir, model_name, "_Misspec_Length_Error_perR_withRegression.png"),
  #        p_length_reg, width = 12, height = 6, dpi = 300)
  
  return(merged_tree)
}

# ==== Load Data ====
cat("Loading data files...\n")
simpleTree <- read.csv("../simpleModel/4-analysis/summary/tree_analysis_results.csv")
realisticTree <- read.csv("../realisticModel/4-analysis/summary/tree_analysis_results.csv")
calibratedTree <- read.csv("../calibratedModel/4-analysis/summary/tree_analysis_results.csv")

simpleParam <- read.csv("../simpleModel/4-analysis/summary/var_counts.csv")
realisticParam <- read.csv("../realisticModel/4-analysis/summary/var_counts.csv")
calibratedParam <- read.csv("../calibratedModel/4-analysis/summary/var_counts.csv")

cat("Data loaded successfully!\nStarting analysis...\n\n")

# ==== Run Analyses ====
simple_results <- analyse_tree("simple", simpleTree, simpleParam)
realistic_results <- analyse_tree("complex", realisticTree, realisticParam)
calibrated_results <- analyse_tree("calibrated", calibratedTree, calibratedParam)