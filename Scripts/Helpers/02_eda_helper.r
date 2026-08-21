create_pce_table1 <- function(data) {
  set.seed(20260821)

  binary_variables <- c(
    "prenatal_nicotine",
    "prenatal_opioid",
    "prenatal_alcohol"
  )

  table_data <- data %>%
    filter(prenatal_cannabis %in% c(0, 1)) %>%
    select(-any_of(c("participant_id", "site"))) %>%
    mutate(
      prenatal_cannabis = factor(
        prenatal_cannabis,
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
      by = prenatal_cannabis,
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
        prenatal_alcohol ~ "Prenatal alcohol exposure"
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
          prenatal_nicotine,
          prenatal_opioid,
          prenatal_alcohol
        ) ~ list(workspace = 2e6),
        child_race ~ list(simulate.p.value = TRUE, B = 100000)
      )
    ) %>%
    gtsummary::bold_p(t = 0.05) %>%
    gtsummary::modify_header(
      label = "**Variable**",
      p.value = "**p-value**"
    ) %>%
    gtsummary::modify_caption(
      "**Table 1. Participant characteristics by prenatal cannabis exposure**"
    ) %>%
    gtsummary::modify_footnote(
      p.value = "ANOVA for continuous variables; Fisher's exact test for categorical variables (Monte Carlo simulation for race). Bold indicates p < 0.05."
    )
}
