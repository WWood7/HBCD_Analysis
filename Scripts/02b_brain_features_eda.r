source("Scripts/00_config.r")


library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)


brain_features <- c(
  "cerebral_wm_vol",
  "cortical_gm_vol"
)

brain_qc_indicators <- c(
  "t2w_qc"
)

treatment_definitions <- c(
  tox_and_report = "prenatal_cannabis",
  tox_only = "prenatal_cannabis_tox"
)

df <- read_csv(
  file.path(preprocessed_dir, "preprocessed_df.csv"),
  show_col_types = FALSE
)


output_dir <- file.path(eda_output_dir, "brain_features")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

safe_wilcox_p <- function(data, value_var) {
  observed <- data %>%
    filter(!is.na(.data[[value_var]]))
  if (!all(c(0, 1) %in% unique(observed$treatment_value))) {
    return(NA_real_)
  }
  group_0 <- observed %>%
    filter(treatment_value == 0) %>%
    pull(all_of(value_var))
  group_1 <- observed %>%
    filter(treatment_value == 1) %>%
    pull(all_of(value_var))
  tryCatch(
    stats::wilcox.test(
      x = group_0,
      y = group_1,
      exact = FALSE
    )$p.value,
    error = function(e) NA_real_
  )
}

brain_summary_list <- list()
brain_test_list <- list()
qc_summary_list <- list()

for (definition_name in names(treatment_definitions)) {
  treatment_var <- treatment_definitions[[definition_name]]
  analysis_data <- df %>%
    filter(.data[[treatment_var]] %in% c(0, 1)) %>%
    mutate(
      treatment_value = as.integer(.data[[treatment_var]]),
      pce_group = factor(
        treatment_value,
        levels = c(0, 1),
        labels = c("PCE = 0", "PCE = 1")
      )
    )

  brain_long <- analysis_data %>%
    select(treatment_value, pce_group, all_of(brain_features)) %>%
    pivot_longer(
      cols = all_of(brain_features),
      names_to = "brain_feature",
      values_to = "value"
    )

  brain_summary_list[[definition_name]] <- brain_long %>%
    group_by(brain_feature, treatment_value, pce_group) %>%
    summarise(
      n_group = n(),
      n_observed = sum(!is.na(value)),
      n_missing = sum(is.na(value)),
      percent_missing = 100 * n_missing / n_group,
      mean = mean(value, na.rm = TRUE),
      sd = stats::sd(value, na.rm = TRUE),
      median = stats::median(value, na.rm = TRUE),
      q1 = stats::quantile(value, 0.25, na.rm = TRUE),
      q3 = stats::quantile(value, 0.75, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(treatment_definition = definition_name, .before = 1)

  brain_test_list[[definition_name]] <- bind_rows(lapply(
    brain_features,
    function(feature) {
      data.frame(
        treatment_definition = definition_name,
        brain_feature = feature,
        wilcoxon_p_value = safe_wilcox_p(
          analysis_data,
          feature
        )
      )
    }
  ))

  brain_density_plot <- brain_long %>%
    filter(!is.na(value)) %>%
    ggplot(
      aes(
        x = value,
        color = pce_group,
        fill = pce_group
      )
    ) +
    geom_density(alpha = 0.25, linewidth = 0.9) +
    facet_wrap(
      ~ brain_feature,
      scales = "free",
      labeller = as_labeller(c(
        cerebral_wm_vol = "Cerebral white matter volume",
        cortical_gm_vol = "Cortical gray matter volume"
      ))
    ) +
    labs(
      title = paste("Brain-feature distributions:", definition_name),
      x = "Volume",
      y = "Density",
      color = NULL,
      fill = NULL
    ) +
    theme_minimal(base_size = 12) +
    theme(plot.title.position = "plot")

  ggsave(
    file.path(
      output_dir,
      paste0(definition_name, "_brain_feature_density.png")
    ),
    brain_density_plot,
    width = 10,
    height = 5.5,
    dpi = 300
  )

  brain_violin_plot <- brain_long %>%
    filter(!is.na(value)) %>%
    ggplot(aes(x = pce_group, y = value, fill = pce_group)) +
    geom_violin(trim = FALSE, alpha = 0.45) +
    geom_boxplot(width = 0.16, outlier.alpha = 0.35) +
    facet_wrap(
      ~ brain_feature,
      scales = "free",
      labeller = as_labeller(c(
        cerebral_wm_vol = "Cerebral white matter volume",
        cortical_gm_vol = "Cortical gray matter volume"
      ))
    ) +
    labs(
      title = paste("Brain features by PCE group:", definition_name),
      x = NULL,
      y = "Volume",
      fill = NULL
    ) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "none",
      plot.title.position = "plot"
    )

  ggsave(
    file.path(
      output_dir,
      paste0(definition_name, "_brain_feature_violin_boxplot.png")
    ),
    brain_violin_plot,
    width = 10,
    height = 5.5,
    dpi = 300
  )

  qc_long <- analysis_data %>%
    select(treatment_value, pce_group, all_of(brain_qc_indicators)) %>%
    pivot_longer(
      cols = all_of(brain_qc_indicators),
      names_to = "qc_indicator",
      values_to = "qc_value"
    )

  qc_summary_list[[definition_name]] <- qc_long %>%
    group_by(qc_indicator, treatment_value, pce_group) %>%
    summarise(
      n_group = n(),
      n_observed = sum(!is.na(qc_value)),
      n_missing = sum(is.na(qc_value)),
      percent_missing = 100 * n_missing / n_group,
      mean = mean(qc_value, na.rm = TRUE),
      sd = stats::sd(qc_value, na.rm = TRUE),
      median = stats::median(qc_value, na.rm = TRUE),
      q1 = stats::quantile(qc_value, 0.25, na.rm = TRUE),
      q3 = stats::quantile(qc_value, 0.75, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(treatment_definition = definition_name, .before = 1)

  qc_distribution <- qc_long %>%
    filter(!is.na(qc_value)) %>%
    count(qc_indicator, pce_group, qc_value, name = "n") %>%
    group_by(qc_indicator, pce_group) %>%
    mutate(percent = 100 * n / sum(n)) %>%
    ungroup()

  qc_plot <- ggplot(
    qc_distribution,
    aes(
      x = factor(qc_value),
      y = percent,
      fill = pce_group
    )
  ) +
    geom_col(position = position_dodge(width = 0.8), width = 0.72) +
    facet_wrap(~ qc_indicator, scales = "free_x") +
    labs(
      title = paste("T2w QC distribution:", definition_name),
      x = "QC value",
      y = "Participants within PCE group (%)",
      fill = NULL
    ) +
    theme_minimal(base_size = 12) +
    theme(plot.title.position = "plot")

  ggsave(
    file.path(
      output_dir,
      paste0(definition_name, "_t2w_qc_distribution.png")
    ),
    qc_plot,
    width = 8,
    height = 5.5,
    dpi = 300
  )
}

brain_feature_summary <- bind_rows(brain_summary_list)
brain_feature_tests <- bind_rows(brain_test_list)
brain_qc_summary <- bind_rows(qc_summary_list)

write_csv(
  brain_feature_summary,
  file.path(output_dir, "brain_feature_distribution_summary.csv")
)
write_csv(
  brain_feature_tests,
  file.path(output_dir, "brain_feature_group_tests.csv")
)
write_csv(
  brain_qc_summary,
  file.path(output_dir, "brain_qc_distribution_summary.csv")
)

message("\n========== Brain-feature distributions by PCE group ==========")
print(as_tibble(brain_feature_summary), n = Inf, width = Inf)
message("\n========== Brain-feature Wilcoxon tests ==========")
print(as_tibble(brain_feature_tests), n = Inf, width = Inf)
message("\n========== Brain QC distributions by PCE group ==========")
print(as_tibble(brain_qc_summary), n = Inf, width = Inf)
message("\nOutputs saved to: ", output_dir)