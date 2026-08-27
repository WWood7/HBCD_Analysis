source("Scripts/00_config.r")

library(dplyr)
library(readr)

preprocessed_data <- read_csv(
  file.path(preprocessed_dir, "preprocessed_df.csv"),
  show_col_types = FALSE
)

birthweight_var <- "birth_weight_lbs"
mediator_vars <- "gestational_age"
covariate_vars <- c(
  "prenatal_nicotine", "prenatal_alcohol", "mother_race",
  "mother_ethnicity", "mother_education", "household_income", "site",
  "mother_age_delivery", "food_insecurity", "mother_employment", 
  "hypertension", "preeclampsia", "oligohydramnios"
)
treatment_definitions <- c(
  self_report = "prenatal_cannabis",
  toxicology = "prenatal_cannabis_tox"
)

required_vars <- unique(c(
  birthweight_var,
  mediator_vars,
  covariate_vars,
  unname(treatment_definitions)
))
missing_vars <- setdiff(required_vars, names(preprocessed_data))
if (length(missing_vars) > 0L) {
  stop("Missing columns: ", paste(missing_vars, collapse = ", "))
}

summarise_birthweight_missingness <- function(data, treatment_name, treatment_var) {
  complete_case_vars <- c(treatment_var, mediator_vars, covariate_vars)

  eligible_data <- data %>%
    filter(.data[[treatment_var]] %in% c(0, 1)) %>%
    filter(complete.cases(pick(all_of(complete_case_vars)))) %>%
    mutate(
      birthweight_missing = is.na(.data[[birthweight_var]]),
      treatment_value = as.character(.data[[treatment_var]])
    )

  overall <- eligible_data %>%
    summarise(
      treatment_value = "Overall",
      n_complete_other_variables = n(),
      n_birthweight_observed = sum(!birthweight_missing),
      n_birthweight_missing = sum(birthweight_missing),
      percent_birthweight_missing = 100 * mean(birthweight_missing)
    )

  by_treatment <- eligible_data %>%
    group_by(treatment_value) %>%
    summarise(
      n_complete_other_variables = n(),
      n_birthweight_observed = sum(!birthweight_missing),
      n_birthweight_missing = sum(birthweight_missing),
      percent_birthweight_missing = 100 * mean(birthweight_missing),
      .groups = "drop"
    )

  bind_rows(overall, by_treatment) %>%
    mutate(
      treatment_definition = treatment_name,
      treatment_variable = treatment_var,
      .before = 1
    )
}

birthweight_missingness <- bind_rows(lapply(
  names(treatment_definitions),
  function(treatment_name) {
    summarise_birthweight_missingness(
      preprocessed_data,
      treatment_name,
      treatment_definitions[[treatment_name]]
    )
  }
))

print(birthweight_missingness, n = Inf)
