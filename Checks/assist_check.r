source("Scripts/00_config.r")

library(dplyr)
library(tidyr)

selfreports <- read.csv(
  file.path(rawdata_dir, "PregnancyAndExposures.csv"),
  check.names = FALSE
)

assist_variables <- c(
  assist_v1 = "pex_bm_assistv1_during__use_003",
  assist_v2 = "pex_bm_assistv2_end__use_003"
)

missing_variables <- setdiff(
  c("participant_id", "session_id", unname(assist_variables)),
  names(selfreports)
)
if (length(missing_variables) > 0L) {
  stop("Missing columns: ", paste(missing_variables, collapse = ", "))
}

assist_data <- selfreports %>%
  filter(session_id %in% c("ses-V01", "ses-V02")) %>%
  select(
    participant_id,
    session_id,
    all_of(unname(assist_variables))
  )

clean_assist_values <- function(x) {
  values <- trimws(as.character(x))
  values[is.na(x) | values == ""] <- NA_character_
  values
}

classify_assist_values <- function(x) {
  values <- clean_assist_values(x)
  observed <- values[!is.na(values)]

  if (length(observed) == 0L) return("Missing")
  if (any(observed == "1")) return("Use (1)")
  if (all(observed == "0")) return("No use (0)")
  if (all(observed %in% c("777", "999"))) {
    return("Unclassifiable (777/999)")
  }
  if (all(observed %in% c("0", "777", "999"))) {
    return("No use + unclassifiable")
  }
  "Other/mixed"
}

# Record-level distributions show whether missingness is concentrated by
# session and distinguish true missing values from 777/999 responses.
record_level_summary <- assist_data %>%
  pivot_longer(
    cols = all_of(unname(assist_variables)),
    names_to = "assist_version",
    values_to = "raw_value"
  ) %>%
  mutate(
    assist_version = recode(
      assist_version,
      !!!setNames(names(assist_variables), assist_variables)
    ),
    response_status = vapply(
      raw_value,
      classify_assist_values,
      character(1)
    )
  ) %>%
  count(
    session_id,
    assist_version,
    response_status,
    name = "n_records"
  ) %>%
  group_by(session_id, assist_version) %>%
  mutate(percent_records = 100 * n_records / sum(n_records)) %>%
  ungroup()

# Match the participant-level reduction used in 01a_prenatal_cannabis.r.
participant_assist <- assist_data %>%
  group_by(participant_id) %>%
  summarise(
    assist_v1_status = classify_assist_values(
      .data[[assist_variables[["assist_v1"]]]]
    ),
    assist_v2_status = classify_assist_values(
      .data[[assist_variables[["assist_v2"]]]]
    ),
    .groups = "drop"
  )

participant_status_summary <- participant_assist %>%
  pivot_longer(
    cols = c(assist_v1_status, assist_v2_status),
    names_to = "assist_version",
    values_to = "response_status"
  ) %>%
  mutate(
    assist_version = recode(
      assist_version,
      assist_v1_status = "assist_v1",
      assist_v2_status = "assist_v2"
    )
  ) %>%
  count(assist_version, response_status, name = "n_participants") %>%
  group_by(assist_version) %>%
  mutate(
    percent_participants = 100 * n_participants / sum(n_participants)
  ) %>%
  ungroup()

participant_availability <- participant_assist %>%
  summarise(
    n_participants = n(),
    both_missing = sum(
      assist_v1_status == "Missing" &
        assist_v2_status == "Missing"
    ),
    v1_only_available = sum(
      assist_v1_status != "Missing" &
        assist_v2_status == "Missing"
    ),
    v2_only_available = sum(
      assist_v1_status == "Missing" &
        assist_v2_status != "Missing"
    ),
    both_available = sum(
      assist_v1_status != "Missing" &
        assist_v2_status != "Missing"
    )
  )

v1_v2_status_cross_tab <- participant_assist %>%
  count(
    assist_v1_status,
    assist_v2_status,
    name = "n_participants"
  ) %>%
  arrange(desc(n_participants))

message("\n========== ASSIST record-level response status ==========")
print(record_level_summary, n = Inf, width = Inf)

message("\n========== ASSIST participant-level response status ==========")
print(participant_status_summary, n = Inf, width = Inf)

message("\n========== ASSIST participant-level availability ==========")
print(participant_availability, width = Inf)

message("\n========== ASSIST V1 by V2 participant status ==========")
print(v1_v2_status_cross_tab, n = Inf, width = Inf)
