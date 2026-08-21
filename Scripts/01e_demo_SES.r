# preprocess all the other variables

source("Scripts/00_config.r")
source("Scripts/Helpers/01_preprocess_helper.r")

# load the libraries
library(dplyr)

# read in the raw tables
demographics <- read.csv(
  file.path(rawdata_dir, "Demographics.csv"),
  check.names = FALSE
)
physicalhealth <- read.csv(
  file.path(rawdata_dir, "PhysicalHealth.csv"),
  check.names = FALSE
)
pregnancyandexposures <- read.csv(
  file.path(rawdata_dir, "PregnancyAndExposures.csv"),
  check.names = FALSE
)
df <- demographics %>%
  full_join(physicalhealth, by = c("participant_id", "session_id")) %>%
  full_join(pregnancyandexposures, by = c("participant_id", "session_id")) %>%
  filter(session_id %in% c("ses-V01", "ses-V02"))





# ------------------------------------------------------------------------------
# Child-related variables
# gestational age at birth: sed_basic_demographics_gestational_age_delivery
# weight at birth: pex_bm_healthv2_inf_001__02 (pounds)
# sex: sed_basic_demographics_sex (0: female, 1: male, 2: unknown)
# child_ethnicity:sed_basic_demographics_child_ethnicity
# child_race: sed_basic_demographics_child_race
# head circumference at scan:
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Parent-related variables
# household income: sed_basic_demographics_rc_mother_income
# mother education: sed_basic_demographics_rc_mother_education
# recruitment site: sed_basic_demographics_recruitment_site
# prenatal other stimulant exposre:
# employment:
# food insecurity:
# insurance type:
# ------------------------------------------------------------------------------

# Variables with names specified above.
demo_ses_variables <- c(
  "sed_basic_demographics_gestational_age_delivery",
  "pex_bm_healthv2_inf_001__02",
  "sed_basic_demographics_sex",
  "sed_basic_demographics_child_ethnicity",
  "sed_basic_demographics_child_race",
  "sed_basic_demographics_rc_mother_income",
  "sed_basic_demographics_rc_mother_education",
  "sed_basic_demographics_recruitment_site"
)

# Reduce V01 and V02 to one record per participant. V01 is checked first, then
# V02; the first non-missing value is retained for each variable independently.
demo_ses_table <- df %>%
  arrange(
    participant_id,
    match(session_id, c("ses-V01", "ses-V02"))
  ) %>%
  group_by(participant_id) %>%
  summarise(
    across(all_of(demo_ses_variables), first_non_missing),
    .groups = "drop"
  )

  demo_ses_table <- demo_ses_table %>%
  apply_variable_spec(demo_ses_spec)

  write.csv(demo_ses_table, file.path(preprocessed_dir, "demo_ses_table.csv"), row.names = FALSE)