create_pce_table1 <- function(
    data,
    group_var,
    definition_name,
    complete_case = FALSE
) {
  set.seed(20260821)

  treatment_variables <- c("prenatal_cannabis", "prenatal_cannabis_tox")
  if (length(group_var) != 1L || !group_var %in% treatment_variables) {
    stop(
      "group_var must be either prenatal_cannabis or prenatal_cannabis_tox."
    )
  }
  if (length(definition_name) != 1L || is.na(definition_name)) {
    stop("definition_name must be one non-missing value.")
  }
  if (length(complete_case) != 1L || is.na(complete_case) ||
      !is.logical(complete_case)) {
    stop("complete_case must be TRUE or FALSE.")
  }
  missing_variables <- setdiff(treatment_variables, names(data))
  if (length(missing_variables) > 0L) {
    stop("Missing columns: ", paste(missing_variables, collapse = ", "))
  }
  other_treatment_var <- setdiff(treatment_variables, group_var)

  binary_variables <- c(
    "prenatal_nicotine",
    "prenatal_opioid",
    "prenatal_alcohol"
  )

  table_data <- data %>%
    filter(.data[[group_var]] %in% c(0, 1)) %>%
    select(-any_of(c("participant_id", "site", other_treatment_var)))

  if (complete_case) {
    table_data <- table_data %>%
      filter(complete.cases(pick(everything())))
  }

  table_data <- table_data %>%
    mutate(
      "{group_var}" := factor(
        .data[[group_var]],
        levels = c(0, 1),
        labels = c("PCE = 0", "PCE = 1")
      ),
      across(
        any_of(binary_variables),
        ~ factor(.x, levels = c(0, 1), labels = c("No", "Yes"))
      )
    )

  table_data %>%
    gtsummary::tbl_summary(
      by = all_of(group_var),
      statistic = list(
        gtsummary::all_continuous() ~ "{mean} ({sd})",
        gtsummary::all_categorical() ~ "{n} ({p}%)"
      ),
      digits = list(
        gtsummary::all_continuous() ~ 2,
        gtsummary::all_categorical() ~ c(0, 1)
      ),
      missing = "ifany",
      missing_text = "Missing",
      label = list(
        gestational_age ~ "Gestational age",
        birth_weight_lbs ~ "Birth weight (lb)",
        child_sex ~ "Child sex",
        child_ethnicity ~ "Child ethnicity",
        child_race ~ "Child race",
        household_income ~ "Household income",
        mother_education ~ "Mother education",
        prenatal_nicotine ~ "Prenatal nicotine exposure",
        prenatal_opioid ~ "Prenatal opioid exposure",
        prenatal_alcohol ~ "Prenatal alcohol exposure",
        mother_race ~ "Mother race",
        mother_ethnicity ~ "Mother ethnicity",
        mother_age_delivery ~ "Mother age at delivery",
        food_insecurity ~ "Food insecurity",
        mother_employment ~ "Mother employment",
        hypertension ~ "Hypertension",
        preeclampsia ~ "Pre-eclampsia",
        oligohydramnios ~ "Oligohydramnios",
        cerebral_wm_vol ~ "Cerebral white matter volume",
        cortical_gm_vol ~ "Cortical gray matter volume",
        t2w_qc ~ "T2w MRI quality control"
      )
    ) %>%
    gtsummary::add_overall(
      last = TRUE,
      col_label = "**Total**"
    ) %>%
    gtsummary::add_p(
      test = list(
        gtsummary::all_continuous() ~ "oneway.test",
        gtsummary::all_categorical() ~ "fisher.test"
      ),
      test.args = list(
        gtsummary::all_continuous() ~ list(var.equal = TRUE),
        c(
          child_sex,
          child_ethnicity,
          household_income,
          mother_education,
          mother_ethnicity,
          prenatal_nicotine,
          prenatal_opioid,
          prenatal_alcohol,
          hypertension,
          preeclampsia,
          oligohydramnios,
          mother_employment,
          food_insecurity
        ) ~ list(workspace = 2e6),
        c(child_race, mother_race) ~ list(simulate.p.value = TRUE, B = 100000)
      )
    ) %>%
    gtsummary::bold_p(t = 0.05) %>%
    gtsummary::modify_header(
      label = "**Variable**",
      p.value = "**p-value**"
    ) %>%
    gtsummary::modify_caption(
      paste0(
        "**Table 1. Participant characteristics by ",
        gsub("_", " ", definition_name),
        if (complete_case) " (complete-case sample)" else "",
        "**"
      )
    ) %>%
    gtsummary::modify_footnote(
      p.value = "ANOVA for continuous variables; Fisher's exact test for categorical variables (Monte Carlo simulation for race). Bold indicates p < 0.05."
    )
}
