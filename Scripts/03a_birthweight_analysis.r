source("Scripts/00_config.r")
source("Scripts/Helpers/03_mediation_analysis_helper.r")

suppressPackageStartupMessages({
  library(sl3)
  library(origami)
  library(dplyr)
  library(readr)
  library(gt)
})

# Number of participant-level outer cross-fitting folds.
num_outer_folds <- 2L
flag_missingness <- FALSE
runtime_preset <- "fast"
set.seed(123)

pipeline_start <- Sys.time()
message("Pipeline started at: ", format(pipeline_start, "%Y-%m-%d %H:%M:%S"))

# Load the participant-level HBCD analysis table.
input_path <- file.path(preprocessed_dir, "preprocessed_df.csv")
df_preprocessed <- readr::read_csv(input_path, show_col_types = FALSE)

# HBCD analysis variables
y_vars <- "birth_weight_lbs"
m_vars <- "gestational_age"
delta_y_var <- "birth_weight_observed"
id_vars <- "participant_id"
treatment_definitions <- c(
  self_report = "prenatal_cannabis",
  toxicology = "prenatal_cannabis_tox"
)


x_vars <- c(
  "prenatal_nicotine", "prenatal_alcohol", "mother_race", "mother_ethnicity", "mother_education",
  "household_income", "site", "mother_age_delivery", "food_insecurity",
  "mother_employment"
)
# x_vars <- c(
#   "prenatal_nicotine_freq", "prenatal_alcohol_freq", "mother_race", "mother_ethnicity", "mother_education",
#   "household_income", "site", "mother_age_delivery", "work_during_pregnancy", "food_insecurity",
#   "mother_employment", "hypertension", "preeclampsia", "oligohydramnios"
# )

# Outcome-observation indicator: 1 when birth weight is recorded, 0 otherwise.
df_preprocessed[[delta_y_var]] <- as.integer(
  !is.na(df_preprocessed[[y_vars]])
)

# Define classification and regression Super Learners.
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
learners <- list(cls_lrnr = cls_lrnr, reg_lrnr = reg_lrnr)

models_output_dir <- file.path(birthweight_output_dir, "Models")
dir.create(models_output_dir, recursive = TRUE, showWarnings = FALSE)

results_list <- lapply(names(treatment_definitions), function(definition_name) {
  treatment_var <- treatment_definitions[[definition_name]]
  analysis_data <- df_preprocessed %>%
    filter(.data[[treatment_var]] %in% c(0, 1))
  if (!all(c(0, 1) %in% unique(analysis_data[[treatment_var]]))) {
    stop(treatment_var, " does not contain both 0 and 1 after subsetting.")
  }
  message(
    "\n========== HBCD birthweight mediation: ",
    treatment_var,
    " 0 vs 1 (N = ",
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

  output_stem <- paste0(
    "birthweight_",
    definition_name
  )
  saveRDS(
    treatment_results,
    file = file.path(models_output_dir, paste0(output_stem, ".rds"))
  )
  readr::write_csv(
    treatment_results,
    file = file.path(models_output_dir, paste0(output_stem, ".csv"))
  )
  treatment_results
})
names(results_list) <- names(treatment_definitions)
results <- dplyr::bind_rows(results_list)

# Create one HTML table for each treatment definition.
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
      estimate,
      se,
      ci_lower = estimate - z_critical * se,
      ci_upper = estimate + z_critical * se,
      p_value = 2 * stats::pnorm(-abs(estimate / se)),
      standardized_estimate,
      standardized_se,
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

tables_output_dir <- file.path(birthweight_output_dir, "Tables")
dir.create(tables_output_dir, recursive = TRUE, showWarnings = FALSE)

for (definition_name in names(results_list)) {
  treatment_results <- results_list[[definition_name]]
  table_data <- make_table_data(treatment_results)

  effect_table <- table_data %>%
    gt::gt(rowname_col = "effect", groupname_col = "outcome") %>%
    gt::tab_header(
      title = gt::md("**HBCD Birthweight Mediation Effects**"),
      subtitle = paste0(
        treatment_definitions[[definition_name]],
        ": 1 versus 0"
      )
    ) %>%
    gt::tab_spanner(
      label = "Estimated effect (lb)",
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
        " (observed birth weight = ",
        treatment_results$n_outcome_observed[1],
        ", missing birth weight = ",
        treatment_results$n_outcome_missing[1],
        ")."
      )
    )

  gt::gtsave(
    effect_table,
    filename = file.path(
      tables_output_dir,
      paste0("birthweight_", definition_name, "_effects.html")
    )
  )
}


pipeline_end <- Sys.time()
message(
  "\nSaved separate treatment results and combined results to ",
  models_output_dir,
  "\nSaved HTML effect tables to ",
  tables_output_dir,
  "\nPipeline finished at: ", format(pipeline_end, "%Y-%m-%d %H:%M:%S"),
  " | Total runtime: ",
  round(as.numeric(difftime(pipeline_end, pipeline_start, units = "mins")), 2),
  " mins"
)