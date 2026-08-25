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
SES <- read.csv(
  file.path(rawdata_dir, "SocialEnvironmentalDeterminants.csv"),
  check.names = FALSE
)

df <- demographics %>%
  full_join(physicalhealth, by = c("participant_id", "session_id")) %>%
  full_join(pregnancyandexposures, by = c("participant_id", "session_id")) %>%
  full_join(SES, by = c("participant_id", "session_id"))





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
# mother race: sed_basic_demographics_screen_mother_race
# mother ethnicity: sed_basic_demographics_screen_mother_ethnicity
# mother age at delivery: sed_basic_demographics_mother_age_delivery
# prenatal other stimulant exposre:
# mother employment: sed_bm_demo_work_001
# food insecurity: sed_cg_foodins_category
# work during pregnancy: sed_cg_employ_001
# insurance type:
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Placental Pathology Proxies
# hypertension: pex_bm_healthv2_preg__compl_001___2
# preeclampsia: pex_bm_healthv2_preg__compl_001___3
# oligohydramnios: pex_bm_healthv2_preg__compl_001___9
# ------------------------------------------------------------------------------

demo_ses_variables <- c(
  "sed_basic_demographics_gestational_age_delivery", # gestational age at birth (weeks)
  "pex_bm_healthv2_inf_001__02", # weight at birth (pounds)
  "sed_basic_demographics_sex", # sex (0: female, 1: male, 2: unknown)
  "sed_basic_demographics_child_ethnicity", # child ethnicity
  "sed_basic_demographics_child_race", # child race
  "sed_basic_demographics_rc_mother_income", # household income
  "sed_basic_demographics_rc_mother_education", # mother education
  "sed_basic_demographics_recruitment_site", # recruitment site
  "sed_basic_demographics_screen_mother_race", # mother race
  "sed_basic_demographics_screen_mother_ethnicity", # mother ethnicity
  "sed_basic_demographics_mother_age_delivery", # mother age at delivery (years)
  "sed_cg_employ_001", # work during pregnancy
  "sed_cg_foodins_category", # food insecurity
  "sed_bm_demo_work_001", # mother employment
  "pex_bm_healthv2_preg__compl_001___2", # hypertension
  "pex_bm_healthv2_preg__compl_001___3", # preeclampsia
  "pex_bm_healthv2_preg__compl_001___9" # oligohydramnios
)

# Reduce V01 and V02 to one record per participant. V02 is checked first, then
# V01; the first non-missing value is retained for each variable independently.
demo_ses_table <- df %>%
  arrange(
    participant_id,
    match(session_id, c("ses-V02", "ses-V01", "ses-V03", "ses-V04", "ses-V05"))
  ) %>%
  group_by(participant_id) %>%
  summarise(
    across(
      all_of(demo_ses_variables),
      ~ {
        variable <- cur_column()
        if (variable %in% c(
          "sed_basic_demographics_sex",
          "sed_basic_demographics_child_ethnicity"
        )) {
          first_preferred_value(.x, values_to_deprioritize = 2)
        } else if (variable == "sed_basic_demographics_child_race") {
          first_preferred_value(.x, values_to_deprioritize = 7)
        } else {
          first_non_missing(.x)
        }
      }
    ),
    .groups = "drop"
  )

  demo_ses_table <- demo_ses_table %>%
  apply_variable_spec(demo_ses_spec)

  write.csv(demo_ses_table, file.path(preprocessed_dir, "demo_ses_table.csv"), row.names = FALSE)