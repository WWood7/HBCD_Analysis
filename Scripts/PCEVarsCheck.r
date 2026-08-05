# read in the data
df <- read.csv(
  "R:/brisk/HBCD/ASSIT_TBFL_Urine_Nail_V1V2.csv",
  check.names = FALSE
)

# Check the availability of nail results
# determine if we can safely drop the toe nail results
library(dplyr)

nail_collection <- df %>%
  transmute(
    participant_id,
    session_id,
    collection_status = bio_bm_biosample_nails_type_collection_status,
    collection_type = bio_bm_biosample_nails_type_collection_nail_type,
    nail_category = case_when(
      collection_status == "" ~ "missing nails",
      collection_status == "Yes" & collection_type == 0 ~ "fingernails",
      collection_status == "Yes" & collection_type == 1 ~ "toenails",
      TRUE ~ NA_character_
    )
  )

# Within each session, show the distribution of nail collection outcomes
nail_proportion_by_session <- nail_collection %>%
  filter(!is.na(nail_category)) %>%
  mutate(
    nail_category = factor(
      nail_category,
      levels = c("missing nails", "toenails", "fingernails")
    )
  ) %>%
  count(session_id, nail_category, .drop = FALSE, name = "n") %>%
  group_by(session_id) %>%
  mutate(proportion = n / sum(n)) %>%
  ungroup()

print(nail_proportion_by_session)

# Aggregate sessions: a participant has fingernails collected if this
# occurred in at least one session
finger_collection_by_participant <- nail_collection %>%
  group_by(participant_id) %>%
  summarise(
    finger_collection = if_else(
      any(nail_category == "fingernails", na.rm = TRUE),
      "has finger nails collected",
      "has no finger nails collected"
    ),
    .groups = "drop"
  )

finger_collection_proportion <- finger_collection_by_participant %>%
  count(finger_collection, name = "n") %>%
  mutate(proportion = n / sum(n))

print(finger_collection_proportion)


# derive a PCE variable
# positive confirmatory test results for urine or nail samples
# Self-report weekly use of cannabis for four weeks or more during pregnancy

tlfb_variables <- sprintf("pex_ch_tlfb_cbd_wk_%02d", 3:9)
nail_thc_variable <-
  "bio_bm_biosample_nails_results_c_delta-9-THC_n_cat"
nail_type_variable <-
  "bio_bm_biosample_nails_type_collection_nail_type"
urine_thc_variable <-
  "bio_bm_biosample_urine_results_bio_c_delta-9-THC_u_cat"
assist_v1_variable <- "pex_bm_assistv1_during__use_003"
assist_v2_variable <- "pex_bm_assistv2_end__use_003"

is_missing_value <- function(x) {
  is.na(x) | trimws(as.character(x)) == ""
}

has_code <- function(x, code) {
  any(trimws(as.character(x)) == as.character(code), na.rm = TRUE)
}

has_positive_count <- function(x) {
  values <- suppressWarnings(as.numeric(trimws(as.character(x))))
  any(values > 0, na.rm = TRUE)
}

all_existing_are <- function(x, codes) {
  observed <- trimws(as.character(x[!is_missing_value(x)]))
  length(observed) > 0 && all(observed %in% as.character(codes))
}

all_missing_or_code <- function(x, codes) {
  values <- trimws(as.character(x))
  all(is_missing_value(x) | values %in% as.character(codes))
}

# Reduce all sessions to one row per participant. Each TLFB week is counted
# once if use was reported for that week in any session.
pce_components <- df %>%
  group_by(participant_id) %>%
  summarise(
    across(all_of(tlfb_variables), has_positive_count),
    nail_positive = has_code(.data[[nail_thc_variable]], 1),
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
    urine_negative = all_existing_are(.data[[urine_thc_variable]], 4),
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

# The explicit missing-biomarker rule is represented by the final
# partial-negative category unless the unclassifiable rule also applies.
pce_proportions <- pce_components %>%
  count(pce_label, .drop = FALSE, name = "n") %>%
  mutate(proportion = n / sum(n))

print(pce_proportions)

# Among more-than-minimal users, identify participants classified solely by
# a positive nail result, then count those whose positive result was from a
# toenail sample.
nail_only_toenail_summary <- pce_components %>%
  filter(pce == 1) %>%
  summarise(
    n_more_than_minimal = n(),
    n_determined_solely_by_nails = sum(
      nail_positive & !urine_positive & tlfb_weeks_used < 4
    ),
    n_determined_solely_by_toenails = sum(
      nail_positive_toenail & !urine_positive & tlfb_weeks_used & !fingernail_positive < 4
    ),
    proportion_determined_solely_by_toenails =
      n_determined_solely_by_toenails / n_more_than_minimal
  )

print(nail_only_toenail_summary)
