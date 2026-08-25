# use only toxicological resutls for PCE definition
# positive confirmatory test results for urine or nail samples

source("Scripts/00_config.r")
source("Scripts/Helpers/01_preprocess_helper.r")

# load the libraries
library(dplyr)
library(tidyr)

# read in the raw tables
biospecimens <- read.csv(
  file.path(rawdata_dir, "BioSpecimens.csv"),
  check.names = FALSE
)
selfreports <- read.csv(
  file.path(rawdata_dir, "PregnancyAndExposures.csv"),
  check.names = FALSE
)

# merge by both participant_id and session_id
# only keep the records from session V1 and V2
df <- biospecimens %>%
  full_join(selfreports, by = c("participant_id", "session_id")) %>%
  filter(session_id %in% c("ses-V01", "ses-V02"))


# list all the needed variables for deriving beyond-threshold prenatal cannabis exposure
tlfb_thc_variables <- sprintf("pex_ch_tlfb_thc_wk_%02d", 3:9)
# nail_thc_variable <-
#   "bio_bm_biosample_nails_results_c_delta-9-THC_n_cat"
# urine_thc_variable <-
#   "bio_bm_biosample_urine_results_bio_c_delta-9-THC_u_cat"
nail_thc_variable <-
  "bio_bm_biosample_nails_results_c_delta-9-THC_n_cat"
urine_thc_variable <-
  "bio_bm_biosample_urine_results_bio_c_delta-9-THC_u_cat"

assist_v1_variable <- "pex_bm_assistv1_during__use_003"
assist_v2_variable <- "pex_bm_assistv2_end__use_003"


# Reduce all sessions to one row per participant. Each TLFB week is counted
# once if use was reported for that week in any session.
pce_components <- df %>%
  group_by(participant_id) %>%
  summarise(
    across(all_of(tlfb_variables), has_positive_count),
    nail_positive = has_code(.data[[nail_thc_variable]], 1),
    urine_positive = has_code(.data[[urine_thc_variable]], 1),
    nail_negative = all_existing_are(.data[[nail_thc_variable]], c(0, 4)),
    urine_negative = all_existing_are(.data[[urine_thc_variable]], c(0, 4)),
    nail_missing = all(is_missing_value(.data[[nail_thc_variable]])),
    urine_missing = all(is_missing_value(.data[[urine_thc_variable]])),
    assist_v1_no_use = all_existing_are(.data[[assist_v1_variable]], 0),
    assist_v2_no_use = all_existing_are(.data[[assist_v2_variable]], 0),
    assist_any_use = has_code(.data[[assist_v1_variable]], 1) |
      has_code(.data[[assist_v2_variable]], 1),
    assist_v1_unclassifiable = all_missing_or_code(
      .data[[assist_v1_variable]], c(777, 999)
    ),
    assist_v2_unclassifiable = all_missing_or_code(
      .data[[assist_v2_variable]], c(777, 999)
    ),
    nail_unclassifiable = all_missing_or_code(
      .data[[nail_thc_variable]], 3
    ),
    urine_unclassifiable = all_missing_or_code(
      .data[[urine_thc_variable]], 3
    ),
    .groups = "drop"
  ) %>%
  mutate(
    tlfb_weeks_used = rowSums(across(all_of(tlfb_variables))),
    pce = case_when(
      tlfb_weeks_used >= 4 | nail_positive | urine_positive ~ 1L,
      (nail_negative | urine_negative) & !nail_positive & !urine_positive &
        assist_v1_no_use & assist_v2_no_use ~ 0L,
      assist_any_use ~ 2L,
      assist_v1_unclassifiable & assist_v2_unclassifiable &
        nail_unclassifiable & urine_unclassifiable ~ 4L,
      TRUE ~ 3L
    ),
    pce_label = factor(
      pce,
      levels = 0:4,
      labels = c(
        "true negative",
        "more than minimal users",
        "below-threshold users",
        "partial negative",
        "unclassifiable"
      )
    ),
    pce_tox = case_when(
        nail_positive | urine_positive ~ 1L,
        (nail_negative | urine_negative) & !nail_positive & !urine_positive ~ 0L,
        TRUE ~ NA_integer_
      )
  ) %>%
  select(-all_of(tlfb_variables))


# now create a new table with just the PCE variable
# unclassifiable participants will be recorded as missing
# more than minimal users will be recorded as 1
# all other participants will be recorded as 0
pce_table <- pce_components %>%
  mutate(
    prenatal_cannabis = pce,
    prenatal_cannabis_tox = case_when(
      pce_tox == 1 ~ 1L,
      pce_tox == 0 ~ 0L,
      TRUE ~ NA_integer_
    )
  ) %>%
  select(participant_id, prenatal_cannabis, prenatal_cannabis_tox)

# write the output to a csv file
write.csv(pce_table, file.path(preprocessed_dir, "prenatal_cannabis.csv"), row.names = FALSE)



