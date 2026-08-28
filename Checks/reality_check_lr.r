source("Scripts/00_config.r")

library(dplyr)
library(readr)

# Use the same outcome, mediator, treatment definitions, and covariates as the
# birthweight mediation analysis, but fit ordinary linear regression models.
outcome_var <- "birth_weight_lbs"
mediator_vars <- "gestational_age"
treatment_definitions <- c(
  tox_and_report = "prenatal_cannabis",
  tox_only = "prenatal_cannabis_tox"
)
covariate_vars <- c(
  "prenatal_nicotine", "prenatal_alcohol", "mother_race",
  "mother_ethnicity", "mother_education", "household_income", "site",
  "mother_age_delivery", "food_insecurity", "mother_employment",
  "hypertension", "preeclampsia", "oligohydramnios"
)

preprocessed_data <- read_csv(
  file.path(preprocessed_dir, "preprocessed_df.csv"),
  show_col_types = FALSE
)

required_vars <- unique(c(
  outcome_var,
  mediator_vars,
  covariate_vars,
  unname(treatment_definitions)
))
missing_vars <- setdiff(required_vars, names(preprocessed_data))
if (length(missing_vars) > 0L) {
  stop("Missing columns: ", paste(missing_vars, collapse = ", "))
}

extract_treatment_result <- function(
    model,
    treatment_var,
    treatment_definition,
    model_name
) {
  coefficient_table <- summary(model)$coefficients
  if (!treatment_var %in% rownames(coefficient_table)) {
    stop("Treatment coefficient is not estimable for ", treatment_definition, ".")
  }

  estimate <- coefficient_table[treatment_var, "Estimate"]
  standard_error <- coefficient_table[treatment_var, "Std. Error"]
  model_summary <- summary(model)

  data.frame(
    treatment_definition = treatment_definition,
    model = model_name,
    n = stats::nobs(model),
    treatment_coefficient = estimate,
    standard_error = standard_error,
    treatment_p_value = coefficient_table[treatment_var, "Pr(>|t|)"],
    r_squared = model_summary$r.squared,
    adjusted_r_squared = model_summary$adj.r.squared
  )
}

linear_models <- list()
linear_regression_results <- lapply(
  names(treatment_definitions),
  function(definition_name) {
    treatment_var <- treatment_definitions[[definition_name]]
    complete_case_vars <- unique(c(
      outcome_var,
      treatment_var,
      mediator_vars,
      covariate_vars
    ))

    # Use one common complete-case sample so the two model estimates are
    # directly comparable within each treatment definition.
    analysis_data <- preprocessed_data %>%
      filter(.data[[treatment_var]] %in% c(0, 1)) %>%
      select(all_of(complete_case_vars)) %>%
      filter(complete.cases(pick(everything()))) %>%
      mutate(
        "{treatment_var}" := as.numeric(.data[[treatment_var]]),
        site = factor(site)
      )

    if (!all(c(0, 1) %in% unique(analysis_data[[treatment_var]]))) {
      stop(
        treatment_var,
        " does not contain both treatment levels after complete-case filtering."
      )
    }

    treatment_only_model <- stats::lm(
      stats::reformulate(treatment_var, response = outcome_var),
      data = analysis_data
    )
    all_variables_model <- stats::lm(
      stats::reformulate(
        c(treatment_var, mediator_vars, covariate_vars),
        response = outcome_var
      ),
      data = analysis_data
    )

    linear_models[[definition_name]] <<- list(
      treatment_only = treatment_only_model,
      all_variables = all_variables_model
    )

    message(
      "\n========== Linear regression: ",
      definition_name,
      " (N = ",
      nrow(analysis_data),
      ") =========="
    )
    message("\nTreatment-only model")
    print(summary(treatment_only_model))
    message("\nAll-variables model")
    print(summary(all_variables_model))

    bind_rows(
      extract_treatment_result(
        treatment_only_model,
        treatment_var,
        definition_name,
        "Treatment only"
      ),
      extract_treatment_result(
        all_variables_model,
        treatment_var,
        definition_name,
        "All variables"
      )
    )
  }
) %>%
  bind_rows()

message("\n========== Treatment coefficient summary ==========")
print(linear_regression_results, row.names = FALSE)
