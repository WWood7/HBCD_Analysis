source("Scripts/00_config.r")

library(dplyr)

brain <- read.csv(
  file.path(rawdata_dir, "Imaging.csv"),
  check.names = FALSE
)

brain_volume_variables <- c(
  left_cerebral_wm =
    "img_bibsnet_space-T2w_desc-aseg_volumes_Left-Cerebral-White-Matter",
  right_cerebral_wm =
    "img_bibsnet_space-T2w_desc-aseg_volumes_Right-Cerebral-White-Matter",
  left_cortical_gm =
    "img_bibsnet_space-T2w_desc-aseg_volumes_Left-Cerebral-Cortex",
  right_cortical_gm =
    "img_bibsnet_space-T2w_desc-aseg_volumes_Right-Cerebral-Cortex"
)

missing_variables <- setdiff(
  c("participant_id", unname(brain_volume_variables)),
  names(brain)
)
if (length(missing_variables) > 0L) {
  stop("Missing columns: ", paste(missing_variables, collapse = ", "))
}

first_non_missing <- function(x) {
  observed <- x[!is.na(x)]
  if (length(observed) == 0L) NA_real_ else observed[[1]]
}

brain_components <- brain %>%
  transmute(
    participant_id,
    left_cerebral_wm = suppressWarnings(as.numeric(as.character(
      .data[[brain_volume_variables[["left_cerebral_wm"]]]]
    ))),
    right_cerebral_wm = suppressWarnings(as.numeric(as.character(
      .data[[brain_volume_variables[["right_cerebral_wm"]]]]
    ))),
    left_cortical_gm = suppressWarnings(as.numeric(as.character(
      .data[[brain_volume_variables[["left_cortical_gm"]]]]
    ))),
    right_cortical_gm = suppressWarnings(as.numeric(as.character(
      .data[[brain_volume_variables[["right_cortical_gm"]]]]
    )))
  ) %>%
  mutate(
    cerebral_wm_vol = left_cerebral_wm + right_cerebral_wm,
    cortical_gm_vol = left_cortical_gm + right_cortical_gm
  )

# Keep one non-missing T2-weighted volume per participant.
brain_table <- brain_components %>%
  group_by(participant_id) %>%
  summarise(
    cerebral_wm_vol = first_non_missing(cerebral_wm_vol),
    cortical_gm_vol = first_non_missing(cortical_gm_vol),
    .groups = "drop"
  )

write.csv(
  brain_table,
  file.path(preprocessed_dir, "brain.csv"),
  row.names = FALSE
)
