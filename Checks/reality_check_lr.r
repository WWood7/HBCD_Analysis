source("Scripts/00_config.r")

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(gt)
})

if (!requireNamespace("car", quietly = TRUE)) {
  stop("Package 'car' is required to calculate VIFs.")
}
if (!requireNamespace("lmerTest", quietly = TRUE)) {
  stop("Package 'lmerTest' is required to fit linear mixed models.")
}

set.seed(20260830)

outcome_var <- "birth_weight_grams"
site_var <- "site"
treatment_definitions <- c(
  tox_and_report = "prenatal_cannabis",
  tox_only = "prenatal_cannabis_tox"
)

prenatal_coexposures <- c(
  "prenatal_nicotine",
  "prenatal_alcohol",
  "prenatal_opioid"
)
baseline_demographics <- c(
  "mother_race",
  "mother_ethnicity",
  "mother_age_delivery"
)
ses_variables <- c(
  "mother_education",
  "household_income",
  "food_insecurity",
  "mother_employment"
)
potential_mediator <- "gestational_age"
pregnancy_conditions <- c(
  "hypertension",
  "preeclampsia",
  "oligohydramnios"
)

# Gestational age is entered first to show its immediate impact on the PCE
# coefficient. Gestational age and pregnancy conditions may change the estimand.
adjustment_sets <- list(
  "PCE only" = character(0),
  "+ Gestational age" = potential_mediator,
  "+ Prenatal co-exposures" = c(
    potential_mediator,
    prenatal_coexposures
  ),
  "+ Baseline demographics" = c(
    potential_mediator,
    prenatal_coexposures,
    baseline_demographics
  ),
  "+ SES" = c(
    potential_mediator,
    prenatal_coexposures,
    baseline_demographics,
    ses_variables
  ),
  "+ Pregnancy conditions" = c(
    potential_mediator,
    prenatal_coexposures,
    baseline_demographics,
    ses_variables,
    pregnancy_conditions
  )
)
adjustment_roles <- c(
  "PCE only" = "Site random-intercept model",
  "+ Gestational age" = "Potential mediator adjustment",
  "+ Prenatal co-exposures" = "Prenatal co-exposure adjustment",
  "+ Baseline demographics" = "Baseline confounder adjustment",
  "+ SES" = "Baseline confounder adjustment",
  "+ Pregnancy conditions" = "Potential post-exposure adjustment"
)

all_adjustment_vars <- unique(unlist(adjustment_sets, use.names = FALSE))
factor_variables <- c(
  "mother_race",
  "mother_ethnicity",
  "site",
  "mother_education",
  "household_income",
  "food_insecurity",
  "mother_employment",
  "hypertension",
  "preeclampsia",
  "oligohydramnios",
  "child_sex"
)
continuous_predictors <- c(
  prenatal_coexposures,
  "mother_age_delivery",
  potential_mediator
)

preprocessed_data <- read_csv(
  file.path(preprocessed_dir, "preprocessed_df.csv"),
  show_col_types = FALSE
)

required_vars <- unique(c(
  outcome_var,
  site_var,
  all_adjustment_vars,
  unname(treatment_definitions)
))
missing_vars <- setdiff(required_vars, names(preprocessed_data))
if (length(missing_vars) > 0L) {
  stop("Missing columns: ", paste(missing_vars, collapse = ", "))
}

output_dir <- file.path(getwd(), "Outputs", "Checks", "reality_check_lr")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

extract_treatment_result <- function(
    model,
    treatment_var,
    treatment_definition,
    model_name,
    model_order,
    adjustment_role
) {
  coefficient_table <- summary(model)$coefficients
  if (!treatment_var %in% rownames(coefficient_table)) {
    stop("Treatment coefficient is not estimable for ", treatment_definition, ".")
  }

  estimate <- coefficient_table[treatment_var, "Estimate"]
  standard_error <- coefficient_table[treatment_var, "Std. Error"]
  degrees_freedom <- if ("df" %in% colnames(coefficient_table)) {
    coefficient_table[treatment_var, "df"]
  } else {
    Inf
  }
  critical_value <- if (is.finite(degrees_freedom)) {
    stats::qt(0.975, df = degrees_freedom)
  } else {
    stats::qnorm(0.975)
  }
  p_value_column <- grep(
    "^Pr\\(",
    colnames(coefficient_table),
    value = TRUE
  )
  if (length(p_value_column) != 1L) {
    stop("Unable to identify the mixed-model p-value column.")
  }

  data.frame(
    treatment_definition = treatment_definition,
    model = model_name,
    model_order = model_order,
    adjustment_role = adjustment_role,
    n = stats::nobs(model),
    pce_estimate = estimate,
    standard_error = standard_error,
    ci_lower = estimate - critical_value * standard_error,
    ci_upper = estimate + critical_value * standard_error,
    p_value = coefficient_table[treatment_var, p_value_column],
    pce_vif = calculate_pce_vif(model, treatment_var),
    aic = stats::AIC(model),
    bic = stats::BIC(model)
  )
}

save_plot <- function(plot, filename, width = 9, height = 6) {
  ggplot2::ggsave(
    filename = file.path(output_dir, filename),
    plot = plot,
    width = width,
    height = height,
    dpi = 300
  )
}

save_lmm_diagnostic_panel <- function(model, filename, title) {
  fitted_values <- stats::fitted(model)
  residual_values <- stats::residuals(model)
  standardized_residuals <- residual_values / stats::sigma(model)
  leverage <- tryCatch(
    stats::hatvalues(model),
    error = function(e) rep(NA_real_, length(residual_values))
  )
  cook_distance <- tryCatch(
    stats::cooks.distance(model),
    error = function(e) rep(NA_real_, length(residual_values))
  )

  grDevices::png(
    filename = file.path(output_dir, filename),
    width = 2400,
    height = 1600,
    res = 240
  )
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit({
    graphics::par(old_par)
    grDevices::dev.off()
  })

  graphics::par(
    mfrow = c(2, 3),
    mar = c(4, 4, 3, 1),
    oma = c(0, 0, 2, 0)
  )
  graphics::plot(
    fitted_values,
    residual_values,
    xlab = "Fitted values",
    ylab = "Residuals",
    main = "Residuals vs fitted"
  )
  graphics::abline(h = 0, lty = 2, col = "grey40")
  graphics::lines(
    stats::lowess(fitted_values, residual_values),
    col = "#C8553D",
    lwd = 2
  )
  stats::qqnorm(
    standardized_residuals,
    main = "Normal Q-Q",
    ylab = "Standardized residuals"
  )
  stats::qqline(standardized_residuals, col = "#C8553D", lwd = 2)
  graphics::plot(
    fitted_values,
    sqrt(abs(standardized_residuals)),
    xlab = "Fitted values",
    ylab = expression(sqrt("|standardized residual|")),
    main = "Scale-location"
  )
  graphics::hist(
    standardized_residuals,
    breaks = 30,
    main = "Standardized residuals",
    xlab = "Standardized residual"
  )
  if (any(is.finite(cook_distance))) {
    graphics::plot(
      cook_distance,
      type = "h",
      xlab = "Observation",
      ylab = "Cook's distance",
      main = "Cook's distance"
    )
  } else {
    graphics::plot.new()
    graphics::title("Cook's distance unavailable")
  }
  if (any(is.finite(leverage))) {
    graphics::plot(
      leverage,
      standardized_residuals,
      xlab = "Leverage",
      ylab = "Standardized residuals",
      main = "Residuals vs leverage"
    )
    graphics::abline(h = 0, lty = 2, col = "grey40")
  } else {
    graphics::plot.new()
    graphics::title("Leverage unavailable")
  }
  graphics::mtext(title, outer = TRUE, line = 0, cex = 1.2)
  invisible(NULL)
}

make_partial_residual_data <- function(model, predictors) {
  model_data <- stats::model.frame(model)
  model_coefficients <- lme4::fixef(model)

  pieces <- lapply(predictors, function(predictor) {
    if (!predictor %in% names(model_data) ||
        !predictor %in% names(model_coefficients) ||
        !is.numeric(model_data[[predictor]]) ||
        !is.finite(model_coefficients[[predictor]])) {
      return(NULL)
    }

    data.frame(
      predictor = predictor,
      predictor_value = model_data[[predictor]],
      partial_residual =
        stats::residuals(model) +
        model_coefficients[[predictor]] * model_data[[predictor]]
    )
  })

  dplyr::bind_rows(pieces)
}

calculate_pce_vif <- function(model, treatment_var) {
  fixed_formula <- lme4::nobars(stats::formula(model))
  fixed_terms <- attr(stats::terms(fixed_formula), "term.labels")
  if (length(fixed_terms) == 1L) return(1)

  vif_output <- car::vif(model)
  if (is.matrix(vif_output)) {
    if (!treatment_var %in% rownames(vif_output)) {
      stop("Treatment term is missing from the car::vif() output.")
    }
    return(unname(vif_output[treatment_var, "GVIF"]))
  }

  if (!treatment_var %in% names(vif_output)) {
    stop("Treatment term is missing from the car::vif() output.")
  }
  unname(vif_output[[treatment_var]])
}

assess_collinearity <- function(model, continuous_vars, treatment_definition) {
  design_matrix <- stats::model.matrix(model)
  design_matrix <- design_matrix[
    ,
    colnames(design_matrix) != "(Intercept)",
    drop = FALSE
  ]
  variable_sd <- apply(design_matrix, 2, stats::sd)
  design_matrix <- design_matrix[
    ,
    is.finite(variable_sd) & variable_sd > 0,
    drop = FALSE
  ]
  scaled_design <- scale(design_matrix)
  condition_number <- tryCatch(
    kappa(scaled_design, exact = TRUE),
    error = function(e) Inf
  )

  model_data <- stats::model.frame(model)
  available_continuous <- intersect(continuous_vars, names(model_data))
  max_continuous_correlation <- NA_real_
  if (length(available_continuous) >= 2L) {
    correlation_matrix <- stats::cor(
      model_data[, available_continuous, drop = FALSE],
      use = "pairwise.complete.obs"
    )
    correlation_values <- abs(correlation_matrix[upper.tri(correlation_matrix)])
    if (length(correlation_values) > 0L) {
      max_continuous_correlation <- max(correlation_values, na.rm = TRUE)
    }
  }

  summary_row <- data.frame(
    treatment_definition = treatment_definition,
    condition_number = condition_number,
    maximum_continuous_correlation = max_continuous_correlation
  )

  summary_row
}

sample_summary_list <- list()
raw_difference_list <- list()
model_result_list <- list()
collinearity_summary_list <- list()
linear_models <- list()

for (definition_name in names(treatment_definitions)) {
  treatment_var <- treatment_definitions[[definition_name]]

  eligible_data <- preprocessed_data %>%
    filter(.data[[treatment_var]] %in% c(0, 1)) %>%
    mutate(
      pce_value = as.integer(.data[[treatment_var]]),
      pce_group = factor(
        pce_value,
        levels = c(0, 1),
        labels = c("PCE = 0", "PCE = 1")
      )
    )

  if (!all(c(0, 1) %in% unique(eligible_data$pce_value))) {
    stop(treatment_var, " does not contain both treatment levels.")
  }

  sample_counts <- eligible_data %>%
    group_by(pce_value, pce_group) %>%
    summarise(
      n_total = n(),
      n_birthweight_observed = sum(!is.na(.data[[outcome_var]])),
      n_birthweight_missing = sum(is.na(.data[[outcome_var]])),
      percent_birthweight_missing =
        100 * n_birthweight_missing / n_total,
      .groups = "drop"
    )

  observed_data <- eligible_data %>%
    filter(!is.na(.data[[outcome_var]]))

  outcome_distribution <- observed_data %>%
    group_by(pce_value, pce_group) %>%
    summarise(
      birthweight_mean = mean(.data[[outcome_var]]),
      birthweight_sd = stats::sd(.data[[outcome_var]]),
      birthweight_median = stats::median(.data[[outcome_var]]),
      birthweight_q1 = stats::quantile(.data[[outcome_var]], 0.25),
      birthweight_q3 = stats::quantile(.data[[outcome_var]], 0.75),
      birthweight_minimum = min(.data[[outcome_var]]),
      birthweight_maximum = max(.data[[outcome_var]]),
      .groups = "drop"
    )

  sample_summary_list[[definition_name]] <- sample_counts %>%
    left_join(
      outcome_distribution,
      by = c("pce_value", "pce_group")
    ) %>%
    mutate(treatment_definition = definition_name, .before = 1)

  raw_group_statistics <- observed_data %>%
    group_by(pce_value) %>%
    summarise(
      n = n(),
      mean = mean(.data[[outcome_var]]),
      variance = stats::var(.data[[outcome_var]]),
      .groups = "drop"
    )
  group_0 <- raw_group_statistics %>% filter(pce_value == 0)
  group_1 <- raw_group_statistics %>% filter(pce_value == 1)
  if (group_0$n < 2L || group_1$n < 2L) {
    stop("At least two observed outcomes are required in each PCE group.")
  }

  raw_difference <- group_1$mean - group_0$mean
  raw_difference_se <- sqrt(
    group_1$variance / group_1$n +
      group_0$variance / group_0$n
  )
  welch_df <- raw_difference_se^4 / (
    (group_1$variance / group_1$n)^2 / (group_1$n - 1) +
      (group_0$variance / group_0$n)^2 / (group_0$n - 1)
  )
  raw_critical_value <- stats::qt(0.975, df = welch_df)
  raw_difference_list[[definition_name]] <- data.frame(
    treatment_definition = definition_name,
    mean_pce_0 = group_0$mean,
    sd_pce_0 = sqrt(group_0$variance),
    mean_pce_1 = group_1$mean,
    sd_pce_1 = sqrt(group_1$variance),
    mean_difference_pce_1_minus_0 = raw_difference,
    standard_error = raw_difference_se,
    ci_lower = raw_difference - raw_critical_value * raw_difference_se,
    ci_upper = raw_difference + raw_critical_value * raw_difference_se,
    p_value = 2 * stats::pt(
      -abs(raw_difference / raw_difference_se),
      df = welch_df
    )
  )

  histogram_plot <- ggplot(
    observed_data,
    aes(
      x = .data[[outcome_var]],
      fill = pce_group,
      color = pce_group
    )
  ) +
    geom_histogram(
      aes(y = after_stat(density)),
      bins = 30,
      position = "identity",
      alpha = 0.25
    ) +
    geom_density(linewidth = 1, adjust = 1) +
    labs(
      title = paste("Birthweight distribution:", definition_name),
      x = "Birthweight (lb)",
      y = "Density",
      fill = NULL,
      color = NULL
    ) +
    theme_minimal(base_size = 12) +
    theme(plot.title.position = "plot")
  save_plot(
    histogram_plot,
    paste0(definition_name, "_birthweight_histogram_density.png")
  )

  violin_plot <- ggplot(
    observed_data,
    aes(x = pce_group, y = .data[[outcome_var]], fill = pce_group)
  ) +
    geom_violin(trim = FALSE, alpha = 0.45) +
    geom_boxplot(width = 0.16, outlier.alpha = 0.4) +
    labs(
      title = paste("Birthweight by PCE group:", definition_name),
      x = NULL,
      y = "Birthweight (lb)",
      fill = NULL
    ) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "none",
      plot.title.position = "plot"
    )
  save_plot(
    violin_plot,
    paste0(definition_name, "_birthweight_violin_boxplot.png")
  )

  complete_case_vars <- unique(c(
    outcome_var,
    treatment_var,
    site_var,
    all_adjustment_vars
  ))
  analysis_data <- eligible_data %>%
    select(all_of(complete_case_vars)) %>%
    filter(complete.cases(pick(everything()))) %>%
    mutate(
      "{treatment_var}" := as.numeric(.data[[treatment_var]]),
      across(any_of(factor_variables), factor)
    )

  if (!all(c(0, 1) %in% unique(analysis_data[[treatment_var]]))) {
    stop(
      treatment_var,
      " does not contain both treatment levels after complete-case filtering."
    )
  }
  if (dplyr::n_distinct(analysis_data[[site_var]]) < 2L) {
    stop("At least two sites are required to fit a random intercept.")
  }

  constant_adjustment_vars <- all_adjustment_vars[
    vapply(
      all_adjustment_vars,
      function(variable) dplyr::n_distinct(analysis_data[[variable]]) < 2L,
      logical(1)
    )
  ]
  if (length(constant_adjustment_vars) > 0L) {
    warning(
      "Dropping constant predictors for ",
      definition_name,
      ": ",
      paste(constant_adjustment_vars, collapse = ", ")
    )
  }

  definition_models <- list()
  definition_results <- lapply(seq_along(adjustment_sets), function(model_index) {
    model_name <- names(adjustment_sets)[[model_index]]
    model_adjustments <- setdiff(
      adjustment_sets[[model_index]],
      constant_adjustment_vars
    )
    fixed_effects <- c(treatment_var, model_adjustments)
    mixed_model_formula <- stats::as.formula(paste0(
      outcome_var,
      " ~ ",
      paste(fixed_effects, collapse = " + "),
      " + (1 | ",
      site_var,
      ")"
    ))
    model <- lmerTest::lmer(
      formula = mixed_model_formula,
      data = analysis_data,
      REML = FALSE
    )
    definition_models[[model_name]] <<- model

    extract_treatment_result(
      model = model,
      treatment_var = treatment_var,
      treatment_definition = definition_name,
      model_name = model_name,
      model_order = model_index,
      adjustment_role = adjustment_roles[[model_name]]
    )
  })
  linear_models[[definition_name]] <- definition_models
  model_result_list[[definition_name]] <- bind_rows(definition_results)

  fully_adjusted_model <- definition_models[[length(definition_models)]]
  save_lmm_diagnostic_panel(
    fully_adjusted_model,
    paste0(definition_name, "_fully_adjusted_diagnostics.png"),
    paste("Fully adjusted model:", definition_name)
  )

  partial_residual_data <- make_partial_residual_data(
    fully_adjusted_model,
    continuous_predictors
  )
  if (nrow(partial_residual_data) > 0L) {
    linearity_plot <- ggplot(
      partial_residual_data,
      aes(x = predictor_value, y = partial_residual)
    ) +
      geom_point(alpha = 0.25, size = 1) +
      geom_smooth(
        method = "lm",
        formula = y ~ x,
        se = FALSE,
        color = "#264653",
        linewidth = 0.8
      ) +
      geom_smooth(
        method = "loess",
        formula = y ~ x,
        se = FALSE,
        color = "#C8553D",
        linewidth = 0.8
      ) +
      facet_wrap(~ predictor, scales = "free") +
      labs(
        title = paste("Partial-residual linearity checks:", definition_name),
        subtitle = "Straight line = fitted linear term; curved line = LOESS",
        x = "Predictor value",
        y = "Partial residual"
      ) +
      theme_minimal(base_size = 11) +
      theme(plot.title.position = "plot")
    save_plot(
      linearity_plot,
      paste0(definition_name, "_partial_residual_linearity.png"),
      width = 11,
      height = 7
    )
  }

  collinearity_summary_list[[definition_name]] <- assess_collinearity(
    fully_adjusted_model,
    continuous_predictors,
    definition_name
  )
}

sample_missingness_results <- bind_rows(sample_summary_list)
raw_outcome_results <- bind_rows(raw_difference_list)
sequential_model_results <- bind_rows(model_result_list)
collinearity_results <- bind_rows(collinearity_summary_list)

coefficient_plot_data <- sequential_model_results %>%
  mutate(
    model = factor(model, levels = names(adjustment_sets)),
    treatment_definition = factor(
      treatment_definition,
      levels = names(treatment_definitions)
    )
  )
coefficient_plot <- ggplot(
  coefficient_plot_data,
  aes(
    x = model,
    y = pce_estimate,
    ymin = ci_lower,
    ymax = ci_upper,
    group = treatment_definition
  )
) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey45") +
  geom_line(color = "#264653", linewidth = 0.7) +
  geom_pointrange(color = "#264653") +
  facet_wrap(~ treatment_definition) +
  labs(
    title = "PCE coefficient across sequential adjustment sets",
    x = "Model specification",
    y = "PCE coefficient (lb)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    plot.title.position = "plot"
  )
save_plot(
  coefficient_plot,
  "pce_coefficient_trajectory.png",
  width = 12,
  height = 7
)

pce_vif_plot <- ggplot(
  coefficient_plot_data,
  aes(
    x = model,
    y = pce_vif,
    group = treatment_definition
  )
) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey45") +
  geom_line(color = "#264653", linewidth = 0.7) +
  geom_point(color = "#264653", size = 2.2) +
  facet_wrap(~ treatment_definition) +
  labs(
    title = "PCE VIF across sequential adjustment sets",
    x = "Model specification",
    y = "PCE VIF"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    plot.title.position = "plot"
  )
save_plot(
  pce_vif_plot,
  "pce_vif_trajectory.png",
  width = 12,
  height = 7
)

write_csv(
  sample_missingness_results,
  file.path(output_dir, "sample_missingness_summary.csv")
)
write_csv(
  raw_outcome_results,
  file.path(output_dir, "raw_outcome_difference.csv")
)
write_csv(
  sequential_model_results,
  file.path(output_dir, "sequential_model_results.csv")
)
write_csv(
  collinearity_results,
  file.path(output_dir, "collinearity_summary.csv")
)
write_csv(
  sequential_model_results %>%
    select(treatment_definition, model, model_order, pce_vif),
  file.path(output_dir, "pce_vif_by_model.csv")
)
old_vif_path <- file.path(output_dir, "vif_results.csv")
if (file.exists(old_vif_path)) file.remove(old_vif_path)

model_table_data <- sequential_model_results %>%
  arrange(treatment_definition, model_order) %>%
  transmute(
    Treatment = treatment_definition,
    Model = model,
    `Adjustment role` = adjustment_role,
    N = n,
    `PCE estimate` = sprintf("%.3f", pce_estimate),
    SE = sprintf("%.3f", standard_error),
    `95% CI` = sprintf("[%.3f, %.3f]", ci_lower, ci_upper),
    p = ifelse(p_value < 0.001, "<0.001", sprintf("%.3f", p_value)),
    `PCE VIF` = sprintf("%.3f", pce_vif)
  )

model_table <- model_table_data %>%
  gt::gt(groupname_col = "Treatment") %>%
  gt::tab_header(
    title = "Sequential birthweight linear mixed models"
  ) %>%
  gt::tab_source_note(
    gt::md(
      paste(
        "All models include a site-level random intercept and are fit by",
        "maximum likelihood.",
        "Gestational age is a potential mediator; pregnancy conditions may",
        "also be post-exposure. Coefficient changes after these stages should",
        "not be interpreted as ordinary baseline confounding adjustment."
      )
    )
  )
gt::gtsave(
  model_table,
  filename = file.path(output_dir, "sequential_model_table.html")
)

message("\n========== Sample and birthweight missingness ==========")
print(as_tibble(sample_missingness_results), n = Inf, width = Inf)
message("\n========== Raw PCE=1 minus PCE=0 birthweight difference ==========")
print(as_tibble(raw_outcome_results), n = Inf, width = Inf)
message("\n========== Sequential PCE coefficient results ==========")
print(
  as_tibble(sequential_model_results) %>%
    select(
      treatment_definition,
      model,
      adjustment_role,
      n,
      pce_estimate,
      standard_error,
      ci_lower,
      ci_upper,
      p_value,
      pce_vif
    ),
  n = Inf,
  width = Inf
)
message("\n========== Collinearity screen ==========")
print(as_tibble(collinearity_results), n = Inf, width = Inf)
message("\nPCE VIF is reported for every sequential model above.")
message("\nDiagnostic outputs saved to: ", output_dir)
