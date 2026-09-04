source("Scripts/00_config.r")

library(dplyr)

brain <- read.csv(
  file.path(rawdata_dir, "Imaging.csv"),
  check.names = FALSE
)

brain_variables <- c(
  left_cerebral_wm =
    "img_bibsnet_space-T2w_desc-aseg_volumes_Left-Cerebral-White-Matter",
  right_cerebral_wm =
    "img_bibsnet_space-T2w_desc-aseg_volumes_Right-Cerebral-White-Matter",
  left_cortical_gm =
    "img_bibsnet_space-T2w_desc-aseg_volumes_Left-Cerebral-Cortex",
  right_cortical_gm =
    "img_bibsnet_space-T2w_desc-aseg_volumes_Right-Cerebral-Cortex",
  t2w_qc = "img_brainswipes_xcpd_hash-0f306a2f+0ef9c88a_T2w_summary_QC",
  t2w_age_adjusted = "img_bibsnet_space-T2w_desc-aseg_volumes_adjusted_age"
)


brain_v02 <- brain %>%
  filter(session_id == "ses-V02")

brain_volume_table <- brain_v02 %>%
  transmute(
    participant_id,
    left_cerebral_wm = suppressWarnings(as.numeric(as.character(
      .data[[brain_variables[["left_cerebral_wm"]]]]
    ))),
    right_cerebral_wm = suppressWarnings(as.numeric(as.character(
      .data[[brain_variables[["right_cerebral_wm"]]]]
    ))),
    left_cortical_gm = suppressWarnings(as.numeric(as.character(
      .data[[brain_variables[["left_cortical_gm"]]]]
    ))),
    right_cortical_gm = suppressWarnings(as.numeric(as.character(
      .data[[brain_variables[["right_cortical_gm"]]]]
    ))),
    t2w_age_adjusted_weeks = suppressWarnings(as.numeric(as.character(
      .data[[brain_variables[["t2w_age_adjusted"]]]]
    )))
  ) %>%
  mutate(
    cerebral_wm_vol = left_cerebral_wm + right_cerebral_wm,
    cortical_gm_vol = left_cortical_gm + right_cortical_gm
  ) %>%
  select(participant_id, cerebral_wm_vol, cortical_gm_vol, t2w_age_adjusted_weeks) %>%
  filter(!is.na(cerebral_wm_vol) | !is.na(cortical_gm_vol) | !is.na(t2w_age_adjusted_weeks)) %>%
  distinct()

brain_qc_table <- brain_v02 %>%
  transmute(
    participant_id,
    t2w_qc = suppressWarnings(as.numeric(as.character(
      .data[[brain_variables[["t2w_qc"]]]]
    )))
  ) %>%
  filter(!is.na(t2w_qc)) %>%
  distinct()

if (anyDuplicated(brain_volume_table$participant_id)) {
  stop("Multiple distinct V02 brain-volume records found for a participant.")
}
if (anyDuplicated(brain_qc_table$participant_id)) {
  stop("Multiple distinct V02 T2w QC records found for a participant.")
}

brain_table <- brain_v02 %>%
  distinct(participant_id) %>%
  left_join(brain_volume_table, by = "participant_id") %>%
  left_join(brain_qc_table, by = "participant_id")

write.csv(
  brain_table,
  file.path(preprocessed_dir, "brain.csv"),
  row.names = FALSE
)
