source("Scripts/00_config.r")
source("Scripts/Helpers/03_birthweight_helper.r")


library(sl3)
library(dplyr)
library(readr)
library(gt)

# Number of participant-level outer cross-fitting folds.
num_outer_folds <- 2L
runtime_preset <- "full"
set.seed(123)

pipeline_start <- Sys.time()
message("Pipeline started at: ", format(pipeline_start, "%Y-%m-%d %H:%M:%S"))

# Load the participant-level HBCD analysis table.
input_path <- file.path(preprocessed_dir, "preprocessed_df.csv")
df_preprocessed <- readr::read_csv(input_path, show_col_types = FALSE)

# Define classification and regression Super Learners.
if (runtime_preset == "fast") {
  cls_learners <- list(
    Lrnr_glmnet$new(nfolds = 5, nlambda = 50),
    Lrnr_ranger$new(probability = TRUE, num.trees = 300, min.node.size = 20)
  )
  reg_learners <- list(
    Lrnr_glmnet$new(nfolds = 5, nlambda = 50),
    Lrnr_ranger$new(num.trees = 300, min.node.size = 20)
  )
} else {
  cls_learners <- list(
    Lrnr_glm$new(),
    Lrnr_ranger$new(probability = TRUE),
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
}

cls_stack <- do.call(Stack$new, cls_learners)
reg_stack <- do.call(Stack$new, reg_learners)
learners <- list(
  cls_lrnr = Lrnr_sl$new(
    learners = cls_stack,
    metalearner = Lrnr_nnls$new(eval_function = loss_loglik_binomial)
  ),
  reg_lrnr = Lrnr_sl$new(
    learners = reg_stack,
    metalearner = Lrnr_nnls$new(eval_function = loss_squared_error)
  )
)

message("\n========== HBCD birthweight mediation: prenatal cannabis 0 vs 1 ==========")
results <- run_birthweight_pipeline(
  df_input = df_preprocessed,
  a_var = A_var,
  x_vars = X_vars,
  m_vars = M_vars,
  y_vars = Y_vars,
  learners = learners,
  num_outer_folds = num_outer_folds
)

models_output_dir <- file.path(birthweight_output_dir, "Models")
dir.create(models_output_dir, recursive = TRUE, showWarnings = FALSE)
out_rds <- file.path(models_output_dir, "birthweight_mediation_results.rds")
out_csv <- file.path(models_output_dir, "birthweight_mediation_results.csv")
saveRDS(results, file = out_rds)
readr::write_csv(results, file = out_csv)





# Create a table of effects and standardized effect sizes.
safe_se <- function(variance) {
  ifelse(is.finite(variance) & variance >= 0, sqrt(variance), NA_real_)
}

effect_table_data <- dplyr::bind_rows(lapply(seq_len(nrow(results)), function(i) {
  effect_estimate <- c(results$ate[i], results$nde[i], results$nie[i])
  effect_se <- safe_se(c(results$ate_var[i], results$nde_var[i], results$nie_var[i]))
  standardized_estimate <- c(results$ate_es[i], results$nde_es[i], results$nie_es[i])
  standardized_se <- safe_se(c(
    results$ate_es_var[i],
    results$nde_es_var[i],
    results$nie_es_var[i]
  ))
  critical_value <- stats::qnorm(0.975)

  data.frame(
    outcome = results$y_var[i],
    effect = c("Total effect (ATE)", "Natural direct effect (NDE)", "Natural indirect effect (NIE)"),
    effect_estimate = effect_estimate,
    effect_se = effect_se,
    effect_ci_lower = effect_estimate - critical_value * effect_se,
    effect_ci_upper = effect_estimate + critical_value * effect_se,
    effect_p = 2 * stats::pnorm(-abs(effect_estimate / effect_se)),
    standardized_estimate = standardized_estimate,
    standardized_se = standardized_se,
    standardized_ci_lower = standardized_estimate - critical_value * standardized_se,
    standardized_ci_upper = standardized_estimate + critical_value * standardized_se,
    standardized_p = 2 * stats::pnorm(-abs(standardized_estimate / standardized_se))
  )
}))

effect_table <- effect_table_data %>%
  gt::gt(rowname_col = "effect", groupname_col = "outcome") %>%
  gt::tab_header(
    title = gt::md("**HBCD Birthweight Mediation Effects**"),
    subtitle = "Prenatal cannabis exposure: 1 versus 0"
  ) %>%
  gt::tab_spanner(
    label = "Effect (lb)",
    columns = c(effect_estimate, effect_se, effect_ci_lower, effect_ci_upper, effect_p)
  ) %>%
  gt::tab_spanner(
    label = "Standardized effect size",
    columns = c(
      standardized_estimate,
      standardized_se,
      standardized_ci_lower,
      standardized_ci_upper,
      standardized_p
    )
  ) %>%
  gt::cols_label(
    effect_estimate = "Estimate",
    effect_se = "SE",
    effect_ci_lower = "95% CI lower",
    effect_ci_upper = "95% CI upper",
    effect_p = "p-value",
    standardized_estimate = "Estimate",
    standardized_se = "SE",
    standardized_ci_lower = "95% CI lower",
    standardized_ci_upper = "95% CI upper",
    standardized_p = "p-value"
  ) %>%
  gt::fmt_number(
    columns = c(
      effect_estimate,
      effect_se,
      effect_ci_lower,
      effect_ci_upper,
      standardized_estimate,
      standardized_se,
      standardized_ci_lower,
      standardized_ci_upper
    ),
    decimals = 3
  ) %>%
  gt::fmt(
    columns = c(effect_p, standardized_p),
    fns = function(x) {
      ifelse(
        is.na(x),
        NA_character_,
        ifelse(x < 0.001, "<0.001", formatC(x, format = "f", digits = 3))
      )
    }
  ) %>%
  gt::tab_source_note(
    source_note = paste0(
      "Two-sided Wald tests with 95% confidence intervals; complete-case N = ",
      results$n_complete[1],
      "."
    )
  )

tables_output_dir <- file.path(birthweight_output_dir, "Tables")
dir.create(tables_output_dir, recursive = TRUE, showWarnings = FALSE)
out_table_html <- file.path(tables_output_dir, "birthweight_mediation_effects.html")
out_table_csv <- file.path(tables_output_dir, "birthweight_mediation_effects.csv")
gt::gtsave(effect_table, filename = out_table_html)
readr::write_csv(effect_table_data, file = out_table_csv)