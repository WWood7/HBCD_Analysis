# ------------------------------------------------------------------------------
# directory configuration
# ------------------------------------------------------------------------------
rawdata_dir <- "R:/brisk/HBCD/UntouchedData"
preprocessed_dir <- "R:/brisk/HBCD/PreprocessedData"
eda_output_dir <- file.path(getwd(), "Outputs", "EDA")
birthweight_output_dir <- file.path(getwd(), "Outputs", "Causal", "Birthweight")


# ------------------------------------------------------------------------------
# demographic and SES variables specification list
# ------------------------------------------------------------------------------
demo_ses_spec <- list(
  sed_basic_demographics_gestational_age_delivery = list(
    new_name = "gestational_age",
    type = "numeric"
  ),
  pex_bm_healthv2_inf_001__02 = list(
    new_name = "birth_weight_lbs",
    type = "numeric"
  ),
  sed_basic_demographics_sex = list(
    new_name = "child_sex",
    type = "factor",
    levels = c(0, 1, 2),
    labels = c("Female", "Male", "Unknown")
  ),
  sed_basic_demographics_child_ethnicity = list(
    new_name = "child_ethnicity",
    type = "factor",
    levels = c(0, 1, 2),
    labels = c("Hispanic", "Non-Hispanic", "Unknown")
  ),
  sed_basic_demographics_child_race = list(
    new_name = "child_race",
    type = "factor",
    levels = c(0, 1, 2, 3, 4, 5, 6, 7),
    labels = c("White", "Black", "American Indian/Alaska Native", "Asian", "Native Hawaiian/Other Pacific Islander", "Two or More Races", "Other", "Unknown")
  ),
  sed_basic_demographics_rc_mother_income = list(
    new_name = "household_income",
    type = "factor",
    levels = c(1, 2, 3, 4, 5, 6, 7),
    labels = c("<50k", "<50k", "<50k", "50-100k", "50-100k", ">100k", ">100k")
  ),
  sed_basic_demographics_rc_mother_education = list(
    new_name = "mother_education",
    type = "factor",
    levels = c(1, 2, 3, 4),
    labels = c("Less than high school", "Some college", "Bachelor", "Graduate or professional degree")
  ),
  sed_basic_demographics_recruitment_site = list(
    new_name = "site",
    type = "factor"
  ),
  sed_basic_demographics_screen_mother_race = list(
    new_name = "mother_race",
    type = "factor",
    levels = c(0, 1, 2, 3, 4, 5, 6, 7),
    labels = c("American Indian/Alaska Native", "Asian", "White", "Black", "American Indian/Alaska Native", "Two or More Races", "Two or More Races", "Other")
  ),
  sed_basic_demographics_screen_mother_ethnicity = list(
    new_name = "mother_ethnicity",
    type = "factor",
    levels = c(0, 1),
    labels = c("Hispanic", "Non-Hispanic")
  ),
  sed_basic_demographics_mother_age_delivery = list(
    new_name = "mother_age_delivery",
    type = "numeric"
  ),
  sed_cg_foodins_category = list(
    new_name = "food_insecurity",
    type = "factor",
    levels = c("Negative", "Positive"),
    labels = c("Negative", "Positive"),
    missing_label = "Negative"
  ),
  sed_bm_demo_work_001 = list(
    new_name = "mother_employment",
    type = "factor",
    levels = c(0, 1, 999, 777),
    labels = c("Unemployed", "Employed", "Don't know/Refused", "Don't know/Refused")
  ),
  pex_bm_healthv2_preg__compl_001___2 = list(
    new_name = "hypertension",
    type = "factor",
    levels = c(0, 1),
    labels = c("No", "Yes")
  ),
  pex_bm_healthv2_preg__compl_001___3 = list(
    new_name = "preeclampsia",
    type = "factor",
    levels = c(0, 1),
    labels = c("No", "Yes")
  ),
  pex_bm_healthv2_preg__compl_001___9 = list(
    new_name = "oligohydramnios",
    type = "factor",
    levels = c(0, 1),
    labels = c("No", "Yes")
  )
)
