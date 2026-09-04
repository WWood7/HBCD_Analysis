source("Scripts/00_config.r")
source("Scripts/Helpers/03_mediation_analysis_helper.r")

suppressPackageStartupMessages({
  library(sl3)
  library(origami)
  library(dplyr)
  library(readr)
  library(gt)
})

# Match the birthweight mediation analysis settings.
num_outer_folds <- 2L
flag_missingness <- TRUE
set.seed(123)

pipeline_start <- Sys.time()
message("Pipeline started at: ", format(pipeline_start, "%Y-%m-%d %H:%M:%S"))

input_path <- file.path(preprocessed_dir, "preprocessed_df.csv")
df_preprocessed <- readr::read_csv(input_path, show_col_types = FALSE)

# Brain mediation analysis variables:
# PCE -> birthweight -> infant brain volume.
y_vars <- c(
  "cerebral_wm_vol",
  "cortical_gm_vol"
)
m_vars <- "birth_weight_grams"
delta_y_var <- "brain_outcomes_observed"
brain_qc_var <- "t2w_qc"
brain_qc_threshold <- 0.5
treatment_definitions <- c(
  tox_and_report = "prenatal_cannabis",
  tox_only = "prenatal_cannabis_tox"
)
x_vars <- c(
  "prenatal_nicotine",
  "prenatal_alcohol",
  "prenatal_opioid",
  "mother_race",
  "mother_ethnicity",
  "mother_education",
  "household_income",
  "site",
  "mother_age_delivery",
  "food_insecurity",
  "mother_employment",
  "child_sex",
  "t2w_age_adjusted_weeks"
)

required_columns <- unique(c(
  y_vars,
  m_vars,
  brain_qc_var,
  unname(treatment_definitions),
  x_vars
))
missing_columns <- setdiff(required_columns, names(df_preprocessed))
if (length(missing_columns) > 0L) {
  stop("Missing columns: ", paste(missing_columns, collapse = ", "))
}

# A brain outcome is eligible only when T2w QC is available and passes the
# prespecified threshold. Complete-case preprocessing subsequently requires
# both requested brain outcomes to be non-missing.
df_preprocessed[[delta_y_var]] <- as.integer(
  !is.na(df_preprocessed[[brain_qc_var]]) &
    df_preprocessed[[brain_qc_var]] >= brain_qc_threshold
)

# Use the same Super Learner library as the birthweight analysis.
cls_learners <- list(
  Lrnr_glm$new(),
  Lrnr_ranger$new(),
  Lrnr_earth$new(),
  Lrnr_mean$new(),
  Lrnr_gbm$new()
)
reg_learners <- list(
  Lrnr_glm$new(),
  Lrnr_ranger$new(),
  Lrnr_earth$new(),
  Lrnr_mean$new(),
  Lrnr_gbm$new()
)

cls_stack <- do.call(Stack$new, cls_learners)
reg_stack <- do.call(Stack$new, reg_learners)

cls_lrnr <- Lrnr_sl$new(
  learners = cls_stack,
  metalearner = Lrnr_nnls$new(eval_function = loss_loglik_binomial)
)
reg_lrnr <- Lrnr_sl$new(
  learners = reg_stack,
  metalearner = Lrnr_nnls$new(eval_function = loss_squared_error)
)
learners <- list(
  cls_lrnr = cls_lrnr,
  reg_lrnr = reg_lrnr
)

brain_output_dir <- file.path(getwd(), "Outputs", "Causal", "Brain")
models_output_dir <- file.path(brain_output_dir, "Models")
tables_output_dir <- file.path(brain_output_dir, "Tables")
dir.create(models_output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tables_output_dir, recursive = TRUE, showWarnings = FALSE)

results_list <- lapply(
  names(treatment_definitions),
  function(definition_name) {
    treatment_var <- treatment_definitions[[definition_name]]
    analysis_data <- df_preprocessed %>%
      filter(.data[[treatment_var]] %in% c(0, 1))

    if (!all(c(0, 1) %in% unique(analysis_data[[treatment_var]]))) {
      stop(treatment_var, " does not contain both 0 and 1 after subsetting.")
    }

    message(
      "\n========== HBCD brain mediation: ",
      treatment_var,
      " 0 vs 1 (N before complete-case filtering = ",
      nrow(analysis_data),
      ") =========="
    )

    treatment_results <- run_mediation_pipeline(
      df_input = analysis_data,
      a_var = treatment_var,
      x_vars = x_vars,
      m_vars = m_vars,
      y_vars = y_vars,
      delta_y_var = delta_y_var,
      learners = learners,
      num_outer_folds = num_outer_folds,
      flag_missingness = flag_missingness
    ) %>%
      mutate(treatment_definition = definition_name)

    output_stem <- paste0("brain_", definition_name)
    saveRDS(
      treatment_results,
      file = file.path(models_output_dir, paste0(output_stem, ".rds"))
    )
    readr::write_csv(
      treatment_results,
      file = file.path(models_output_dir, paste0(output_stem, ".csv"))
    )
    treatment_results
  }
)
names(results_list) <- names(treatment_definitions)
results <- dplyr::bind_rows(results_list)

saveRDS(
  results,
  file = file.path(models_output_dir, "brain_combined.rds")
)
readr::write_csv(
  results,
  file = file.path(models_output_dir, "brain_combined.csv")
)

safe_se <- function(variance) {
  ifelse(is.finite(variance) & variance >= 0, sqrt(variance), NA_real_)
}

make_table_data <- function(treatment_results) {
  dplyr::bind_rows(lapply(seq_len(nrow(treatment_results)), function(i) {
    estimate <- c(
      treatment_results$ate[i],
      treatment_results$nde[i],
      treatment_results$nie[i]
    )
    se <- safe_se(c(
      treatment_results$ate_var[i],
      treatment_results$nde_var[i],
      treatment_results$nie_var[i]
    ))
    standardized_estimate <- c(
      treatment_results$ate_es[i],
      treatment_results$nde_es[i],
      treatment_results$nie_es[i]
    )
    standardized_se <- safe_se(c(
      treatment_results$ate_es_var[i],
      treatment_results$nde_es_var[i],
      treatment_results$nie_es_var[i]
    ))
    z_critical <- stats::qnorm(0.975)

    data.frame(
      outcome = treatment_results$y_var[i],
      effect = c(
        "Total effect (ATE)",
        "Natural direct effect (NDE)",
        "Natural indirect effect (NIE)"
      ),
      estimate = estimate,
      se = se,
      ci_lower = estimate - z_critical * se,
      ci_upper = estimate + z_critical * se,
      p_value = 2 * stats::pnorm(-abs(estimate / se)),
      standardized_estimate = standardized_estimate,
      standardized_se = standardized_se,
      standardized_ci_lower =
        standardized_estimate - z_critical * standardized_se,
      standardized_ci_upper =
        standardized_estimate + z_critical * standardized_se,
      standardized_p_value =
        2 * stats::pnorm(-abs(standardized_estimate / standardized_se))
    )
  }))
}

format_p_value <- function(x) {
  ifelse(
    is.na(x),
    NA_character_,
    ifelse(x < 0.001, "<0.001", formatC(x, format = "f", digits = 3))
  )
}

for (definition_name in names(results_list)) {
  treatment_results <- results_list[[definition_name]]
  table_data <- make_table_data(treatment_results)

  effect_table <- table_data %>%
    gt::gt(rowname_col = "effect", groupname_col = "outcome") %>%
    gt::tab_header(
      title = gt::md("**HBCD Brain-Volume Mediation Effects**"),
      subtitle = paste0(
        treatment_definitions[[definition_name]],
        ": 1 versus 0; mediator = birthweight"
      )
    ) %>%
    gt::tab_spanner(
      label = "Estimated effect (volume units)",
      columns = c(estimate, se, ci_lower, ci_upper, p_value)
    ) %>%
    gt::tab_spanner(
      label = "Standardized effect",
      columns = c(
        standardized_estimate,
        standardized_se,
        standardized_ci_lower,
        standardized_ci_upper,
        standardized_p_value
      )
    ) %>%
    gt::cols_label(
      estimate = "Estimate",
      se = "SE",
      ci_lower = "95% CI lower",
      ci_upper = "95% CI upper",
      p_value = "p-value",
      standardized_estimate = "Estimate",
      standardized_se = "SE",
      standardized_ci_lower = "95% CI lower",
      standardized_ci_upper = "95% CI upper",
      standardized_p_value = "p-value"
    ) %>%
    gt::fmt_number(
      columns = c(
        estimate,
        se,
        ci_lower,
        ci_upper,
        standardized_estimate,
        standardized_se,
        standardized_ci_lower,
        standardized_ci_upper
      ),
      decimals = 3
    ) %>%
    gt::fmt(
      columns = c(p_value, standardized_p_value),
      fns = format_p_value
    ) %>%
    gt::tab_source_note(
      source_note = paste0(
        "Two-sided Wald tests and 95% confidence intervals; N = ",
        treatment_results$n_analysis[1],
        ". Complete-case analysis requires T2w QC >= ",
        brain_qc_threshold,
        ", birthweight, and both brain outcomes to be observed."
      )
    )

  gt::gtsave(
    effect_table,
    filename = file.path(
      tables_output_dir,
      paste0("brain_", definition_name, "_effects.html")
    )
  )
}

pipeline_end <- Sys.time()
message(
  "\nSaved separate and combined treatment results to ",
  models_output_dir,
  "\nSaved HTML effect tables to ",
  tables_output_dir,
  "\nPipeline finished at: ",
  format(pipeline_end, "%Y-%m-%d %H:%M:%S"),
  " | Total runtime: ",
  round(
    as.numeric(difftime(pipeline_end, pipeline_start, units = "mins")),
    2
  ),
  " mins"
)
