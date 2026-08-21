# ------------------------------------------------------------------------------
# directory configuration
# ------------------------------------------------------------------------------
rawdata_dir <- "R:/brisk/HBCD/UntouchedData"
preprocessed_dir <- "R:/brisk/HBCD/PreprocessedData"


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
  )
)
