# a large number of participants have unknown sex
# it turns out that many participants with unknown sex recorded at ses-V01
# have no subsequent records
# it seems to be the same reason for unkonwn or missing race and GA too

source("Scripts/00_config.r")
source("Scripts/Helpers/01_preprocess_helper.r")

# load the libraries
library(dplyr)
library(tidyr)

# read in the demographics table
demographics <- read.csv(
  file.path(rawdata_dir, "Demographics.csv"),
  check.names = FALSE
)

# count the number of participants with unknown sex in session 1 (before birth)
# and then check if they have any records in the subsequent sessions
unknown_sex_v01 <- demographics %>%
  filter(
    session_id == "ses-V01",
    trimws(as.character(sed_basic_demographics_sex)) == "2"
  ) %>%
  distinct(participant_id)

# Summarize all post-V01 demographic records for these participants.
followup_by_participant <- unknown_sex_v01 %>%
  left_join(
    demographics %>%
      filter(session_id != "ses-V01") %>%
      semi_join(unknown_sex_v01, by = "participant_id") %>%
      group_by(participant_id) %>%
      summarise(
        n_subsequent_sessions = n_distinct(session_id),
        subsequent_sessions = paste(sort(unique(session_id)), collapse = ", "),
        subsequent_sex_values = paste(
          sort(unique(na.omit(as.character(sed_basic_demographics_child_race)))),
          collapse = ", "
        ),
        .groups = "drop"
      ),
    by = "participant_id"
  ) %>%
  mutate(
    has_subsequent_records = !is.na(n_subsequent_sessions),
    n_subsequent_sessions = coalesce(n_subsequent_sessions, 0L)
  )

followup_summary <- followup_by_participant %>%
  summarise(
    n_unknown_sex_v01 = n(),
    n_with_subsequent_records = sum(has_subsequent_records),
    n_without_subsequent_records = sum(!has_subsequent_records),
    percent_with_subsequent_records = round(
      100 * mean(has_subsequent_records),
      1
    )
  )

subsequent_session_summary <- demographics %>%
  semi_join(unknown_sex_v01, by = "participant_id") %>%
  filter(session_id != "ses-V01") %>%
  distinct(participant_id, session_id) %>%
  count(session_id, name = "n_participants") %>%
  arrange(session_id)

print(followup_summary)
print(subsequent_session_summary)







unknown_race_v01 <- demographics %>%
  filter(
    session_id == "ses-V01",
    trimws(as.character(sed_basic_demographics_child_race)) == "7"
  ) %>%
  distinct(participant_id)

# Summarize all post-V01 demographic records for these participants.
followup_by_participant <- unknown_race_v01 %>%
  left_join(
    demographics %>%
      filter(session_id != "ses-V01") %>%
      semi_join(unknown_race_v01, by = "participant_id") %>%
      group_by(participant_id) %>%
      summarise(
        n_subsequent_sessions = n_distinct(session_id),
        subsequent_sessions = paste(sort(unique(session_id)), collapse = ", "),
        subsequent_sex_values = paste(
          sort(unique(na.omit(as.character(sed_basic_demographics_child_race)))),
          collapse = ", "
        ),
        .groups = "drop"
      ),
    by = "participant_id"
  ) %>%
  mutate(
    has_subsequent_records = !is.na(n_subsequent_sessions),
    n_subsequent_sessions = coalesce(n_subsequent_sessions, 0L)
  )

followup_summary <- followup_by_participant %>%
  summarise(
    n_unknown_race_v01 = n(),
    n_with_subsequent_records = sum(has_subsequent_records),
    n_without_subsequent_records = sum(!has_subsequent_records),
    percent_with_subsequent_records = round(
      100 * mean(has_subsequent_records),
      1
    )
  )

subsequent_session_summary <- demographics %>%
  semi_join(unknown_race_v01, by = "participant_id") %>%
  filter(session_id != "ses-V01") %>%
  distinct(participant_id, session_id) %>%
  count(session_id, name = "n_participants") %>%
  arrange(session_id)

print(followup_summary)
print(subsequent_session_summary)