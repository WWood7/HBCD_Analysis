# Self-report use equal to 7 or more standard drinks per week for two or more weeks during pregnancy
# or
# Self-report use equal to three or more standard drinks per occasion on two or more occasions during pregnancy
# or
# Newborn diagnosed with Fetal Alcohol Syndrome
# or
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
tlfb_variables <- sprintf("pex_ch_tlfb_alc_wk_%02d", 3:9)
tlfb_flag_variable <- "pex_ch_tlfb_self_report_alcohol"
nail_alcohol_variable <-
  "bio_bm_biosample_nails_results_c_ethanol_n"
urine_alcohol_variable <-
  "bio_bm_biosample_urine_results_bio_c_ethanol_u"
infant_fas_variable <-
  "pex_bm_healthv2_inf_007___5"
assist_v1_variable <- "pex_bm_assistv1_during__use_002"
assist_v2_variable <- "pex_bm_assistv2_end__use_002"


# Reduce all sessions to one row per participant. Each TLFB week is counted
# once if use was reported for that week in any session.
pae_components <- df %>%
  group_by(participant_id) %>%
  summarise(
    across(
      all_of(tlfb_variables),
      ~ has_count_at_least(.x, threshold = 7)
    ),
    infant_positive = has_code(.data[[infant_fas_variable]], 1),
    tlfb_flag_positive = has_code(.data[[tlfb_flag_variable]], 1),
    nail_positive = has_code(.data[[nail_alcohol_variable]], 1),
    urine_positive = has_code(.data[[urine_alcohol_variable]], 1),
    nail_negative = all_existing_are(.data[[nail_alcohol_variable]], c(0, 4)),
    urine_negative = all_existing_are(.data[[urine_alcohol_variable]], 0),
    assist_v1_unclassifiable = all_missing_or_code(
      .data[[assist_v1_variable]], c(777, 999)
    ),
    assist_v2_unclassifiable = all_missing_or_code(
      .data[[assist_v2_variable]], c(777, 999)
    ),
    nail_unclassifiable = all_missing_or_code(
      .data[[nail_alcohol_variable]], 3
    ),
    urine_unclassifiable = all_missing_or_code(
      .data[[urine_alcohol_variable]], 3
    ),
    infant_unclassifiable = all_missing_or_code(
      .data[[infant_fas_variable]], NA
    ),
    .groups = "drop"
  ) %>%
  mutate(
    tlfb_weeks_used = rowSums(across(all_of(tlfb_variables))),
    prenatal_alcohol = case_when(
      tlfb_flag_positive | nail_positive | urine_positive | infant_positive ~ 1L,
      assist_v1_unclassifiable & assist_v2_unclassifiable &
        nail_unclassifiable & urine_unclassifiable & infant_unclassifiable ~ NA_integer_,
      TRUE ~ 0L
    )
  ) %>%
  select(-all_of(tlfb_variables))

  # now create a new table with just the PAE variable
  pae_table <- pae_components %>%
    select(participant_id, prenatal_alcohol)

  # write the output to a csv file
  write.csv(pae_table, file.path(preprocessed_dir, "prenatal_alcohol.csv"), row.names = FALSE)



