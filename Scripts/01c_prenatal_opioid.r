# positive toxicology test results
# Self-report use of prescribed or illicit opioids for two weeks or more during pregnancy
# Newborn diagnosed with Neonatal Opioid Withdrawal Syndrome
# Positive toxicology for an opioid in research collected biospecimen

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
physicalhealth <- read.csv(
  file.path(rawdata_dir, "PhysicalHealth.csv"),
  check.names = FALSE
)

# merge by both participant_id and session_id
# only keep the records from session V1 and V2
df <- biospecimens %>%
  full_join(selfreports, by = c("participant_id", "session_id")) %>%
  full_join(physicalhealth, by = c("participant_id", "session_id")) %>%
  filter(session_id %in% c("ses-V01", "ses-V02"))


# list all the needed variables for deriving beyond-threshold prenatal cannabis exposure
tlfb_variables <- sprintf("pex_ch_tlfb_opd_wk_%02d", 3:9)
nail_opioid_variable <-
  "bio_bm_biosample_nails_results_c_any_opioid_n"
urine_opioid_variable <-
  "bio_bm_biosample_urine_results_bio_c_any_opioid_u"
infant_nows_variable <-
  "pex_bm_healthv2_inf_007___1"
assist_v1_variable_prescribed <- "pex_bm_assistv1_during__use_006"
assist_v2_variable_prescribed <- "pex_bm_assistv2_end__use_006"
assist_v1_variable_illicit <- "pex_bm_assistv1_during__use_007"
assist_v2_variable_illicit <- "pex_bm_assistv2_end__use_007"


# Reduce all sessions to one row per participant. Each TLFB week is counted
# once if use was reported for that week in any session.
poe_components <- df %>%
  group_by(participant_id) %>%
  summarise(
    across(all_of(tlfb_variables), has_positive_count),
    infant_positive = has_code(.data[[infant_nows_variable]], 1),
    nail_positive = has_code(.data[[nail_opioid_variable]], 1),
    urine_positive = has_code(.data[[urine_opioid_variable]], 1),
    nail_negative = all_existing_are(.data[[nail_opioid_variable]], c(0, 4)),
    urine_negative = all_existing_are(.data[[urine_opioid_variable]], 0),
    assist_v1_no_use_prescribed = all_existing_are(.data[[assist_v1_variable_prescribed]], 0),
    assist_v2_no_use_prescribed = all_existing_are(.data[[assist_v2_variable_prescribed]], 0),
    assist_v1_no_use_illicit = all_existing_are(.data[[assist_v1_variable_illicit]], 0),
    assist_v2_no_use_illicit = all_existing_are(.data[[assist_v2_variable_illicit]], 0),
    assist_any_use_prescribed = has_code(.data[[assist_v1_variable_prescribed]], 1) |
      has_code(.data[[assist_v2_variable_prescribed]], 1),
    assist_any_use_illicit = has_code(.data[[assist_v1_variable_illicit]], 1) |
      has_code(.data[[assist_v2_variable_illicit]], 1),
    assist_v1_unclassifiable_prescribed = all_missing_or_code(
      .data[[assist_v1_variable_prescribed]], c(777, 999)
    ),
    assist_v2_unclassifiable_prescribed = all_missing_or_code(
      .data[[assist_v2_variable_prescribed]], c(777, 999)
    ),
    assist_v1_unclassifiable_illicit = all_missing_or_code(
      .data[[assist_v1_variable_illicit]], c(777, 999)
    ),
    assist_v2_unclassifiable_illicit = all_missing_or_code(
      .data[[assist_v2_variable_illicit]], c(777, 999)
    ),
    nail_unclassifiable = all_missing_or_code(
      .data[[nail_opioid_variable]], 3
    ),
    urine_unclassifiable = all_missing_or_code(
      .data[[urine_opioid_variable]], 3
    ),
    infant_unclassifiable = all_missing_or_code(
      .data[[infant_nows_variable]], NA
    ),
    .groups = "drop"
  ) %>%
  mutate(
    tlfb_weeks_used = rowSums(across(all_of(tlfb_variables))),
    prenatal_opioid = case_when(
      tlfb_weeks_used >= 2 | nail_positive | urine_positive | infant_positive ~ 1L,
      assist_v1_unclassifiable_prescribed & assist_v2_unclassifiable_prescribed &
        assist_v1_unclassifiable_illicit & assist_v2_unclassifiable_illicit &
        nail_unclassifiable & urine_unclassifiable & infant_unclassifiable ~ NA_integer_,
      TRUE ~ 0L
    )
  ) %>%
  select(-all_of(tlfb_variables))

  # now create a new table with just the PNE variable
  poe_table <- poe_components %>%
    select(participant_id, prenatal_opioid)

  # write the output to a csv file
  write.csv(poe_table, file.path(preprocessed_dir, "prenatal_opioid.csv"), row.names = FALSE)



