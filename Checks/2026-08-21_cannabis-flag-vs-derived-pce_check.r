# a script to check the discrepancy between the source flag and the derived variable
# 8-21-2026 conclusion: the flag is not a good indicator of cannabis exposure
# there are participants with positive with absolute 0 positive indicators

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
demographics <- read.csv(
  file.path(rawdata_dir, "Demographics.csv"),
  check.names = FALSE
)

# merge by both participant_id and session_id
# only keep the records from session V1 and V2
df <- biospecimens %>%
  full_join(selfreports, by = c("participant_id", "session_id")) %>%
  full_join(demographics, by = c("participant_id", "session_id")) %>%
  filter(session_id %in% c("ses-V01", "ses-V02"))


# list all the needed variables for deriving beyond-threshold prenatal cannabis exposure
tlfb_thc_variables <- sprintf("pex_ch_tlfb_thc_wk_%02d", 3:9)
tlfb_cannabinoid_variables <- sprintf("pex_ch_tlfb_cbd_wk_%02d", 3:9)
tlfb_cannibidiol_variables <- sprintf("pex_ch_tlfb_oth_cbd_wk_%02d", 3:9)
tlfb_variables <- sprintf("tlfb_combined_wk_%02d", 3:9)
# Add the cannabinoid and other-cannabinoid values within each week. If both
# source values are missing, keep the combined value missing.
for (i in seq_along(tlfb_variables)) {
  cannabinoid_value <- suppressWarnings(as.numeric(
    trimws(as.character(df[[tlfb_cannabinoid_variables[[i]]]]))
  ))
  other_cannabinoid_value <- suppressWarnings(as.numeric(
    trimws(as.character(df[[tlfb_cannibidiol_variables[[i]]]]))
  ))
  thc_value <- suppressWarnings(as.numeric(
    trimws(as.character(df[[tlfb_thc_variables[[i]]]]))
  ))
  combined_value <- rowSums(
    cbind(cannabinoid_value, other_cannabinoid_value, thc_value),
    na.rm = TRUE
  )
  combined_value[is.na(cannabinoid_value) & is.na(other_cannabinoid_value)] <-
    NA_real_
  df[[tlfb_variables[[i]]]] <- combined_value
}
nail_thc_variable <-
  "bio_bm_biosample_nails_results_c_delta-9-THC_n_cat"
nail_type_variable <-
  "bio_bm_biosample_nails_type_collection_nail_type"
urine_thc_variable <-
  "bio_bm_biosample_urine_results_bio_c_delta-9-THC_u_cat"
assist_v1_variable <- "pex_bm_assistv1_during__use_003"
assist_v2_variable <- "pex_bm_assistv2_end__use_003"
flag_variable <- "par_visit_data_su_flag_cannabis"
flag_tlfb_variable <- "par_visit_data_su_flag_tlfb_bm_cannabis"
tlfb_flag_variable <- "pex_ch_tlfb_self_report_cannabis"


# Reduce all sessions to one row per participant. Each TLFB week is counted
# once if use was reported for that week in any session.
pce_components <- df %>%
  group_by(participant_id) %>%
  summarise(
    across(all_of(tlfb_variables), has_positive_count),
    flag_positive = has_code(.data[[flag_variable]], "Yes"),
    flag_tlfb_positive = has_code(.data[[flag_tlfb_variable]], "Yes"),
    nail_positive = has_code(.data[[nail_thc_variable]], 1),
    nail_screen_positive = has_code(.data[[nail_thc_variable]], c(1, 0)),
    urine_screen_positive = has_code(.data[[urine_thc_variable]], c(1, 0)),
    tlfb_flag_positive = has_code(.data[[tlfb_flag_variable]], 1),
    nail_positive_toenail = any(
      trimws(as.character(.data[[nail_thc_variable]])) == "1" &
        trimws(as.character(.data[[nail_type_variable]])) == "1",
      na.rm = TRUE
    ),
    fingernail_positive = any(
      trimws(as.character(.data[[nail_thc_variable]])) == "1" &
      trimws(as.character(.data[[nail_type_variable]])) == "2",
      na.rm = TRUE
    ),
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
    )
  ) %>%
  select(-all_of(tlfb_variables))









# Review participants whose source flag indicates cannabis exposure but whose
# derived PCE classification is true negative.
flag_positive_pce_zero <- pce_components %>%
  filter(flag_positive, pce != 1L) %>%
  select(
    participant_id,
    nail_positive,
    urine_positive,
    tlfb_weeks_used,
    flag_tlfb_positive,
    nail_screen_positive,
    urine_screen_positive,
    tlfb_flag_positive
  )

print(flag_positive_pce_zero)

flag_negative_pce_one <- pce_components %>%
  filter(!flag_positive, pce == 1L) %>%
  select(
    participant_id,
    nail_positive,
    urine_positive,
    tlfb_weeks_used
  )
print(flag_negative_pce_one)

# Inspect weekly TLFB cannabis responses for participants with a positive TLFB
# flag. Sessions are ordered V01 then V02, and each weekly variable uses its
# first non-missing value across those sessions.
first_non_missing <- function(x) {
  observed_positions <- which(!is_missing_value(x))
  if (length(observed_positions) == 0) {
    return(x[NA_integer_][1])
  }
  x[[observed_positions[[1]]]]
}

tlfb_flagged_records <- df %>%
  arrange(
    participant_id,
    match(session_id, c("ses-V01", "ses-V02"))
  ) %>%
  group_by(participant_id) %>%
  summarise(
    !!flag_tlfb_variable := if_else(
      has_code(.data[[flag_tlfb_variable]], "Yes"),
      "Yes",
      "No"
    ),
    across(all_of(tlfb_variables), first_non_missing),
    .groups = "drop"
  ) %>%
  filter(.data[[flag_tlfb_variable]] == "Yes") %>%
  mutate(
    tlfb_sum = rowSums(
      across(all_of(tlfb_variables)),
      na.rm = TRUE
    )
  ) %>%
  tibble::as_tibble()

print(tlfb_flagged_records, n = Inf, width = Inf)

n_tlfb_flagged_sum_zero <- sum(tlfb_flagged_records$tlfb_sum == 0)
print(n_tlfb_flagged_sum_zero)

# Compare each participant's non-missing TLFB values between V01 and V02.
# Repeated values created by joins are reduced to their distinct value(s)
# within each participant, session, and TLFB variable.
tlfb_observed_by_session <- df %>%
  select(participant_id, session_id, all_of(tlfb_variables)) %>%
  distinct() %>%
  pivot_longer(
    cols = all_of(tlfb_variables),
    names_to = "tlfb_variable",
    values_to = "tlfb_value"
  ) %>%
  filter(!is_missing_value(tlfb_value)) %>%
  mutate(tlfb_value = trimws(as.character(tlfb_value))) %>%
  group_by(participant_id, session_id, tlfb_variable) %>%
  summarise(
    tlfb_value = paste(sort(unique(tlfb_value)), collapse = " | "),
    .groups = "drop"
  )

tlfb_session_comparison <- tlfb_observed_by_session %>%
  filter(session_id == "ses-V01") %>%
  select(participant_id, tlfb_variable, value_v01 = tlfb_value) %>%
  inner_join(
    tlfb_observed_by_session %>%
      filter(session_id == "ses-V02") %>%
      select(participant_id, tlfb_variable, value_v02 = tlfb_value),
    by = c("participant_id", "tlfb_variable")
  ) %>%
  mutate(values_differ = value_v01 != value_v02)

# Cohort summary among participants with non-missing values in both sessions.
tlfb_difference_summary <- tlfb_session_comparison %>%
  group_by(tlfb_variable) %>%
  summarise(
    participants_compared = n(),
    participants_different = sum(values_differ),
    percent_different = 100 * mean(values_differ),
    .groups = "drop"
  )

print(tlfb_difference_summary, n = Inf, width = Inf)