source("Scripts/00_config.r")

library(dplyr)
library(ggplot2)
library(readr)

# Use the participant-level birthweight variable after pounds and ounces have
# been combined in 01e_demo_SES.r.
preprocessed_data <- read_csv(
  file.path(preprocessed_dir, "preprocessed_df.csv"),
  show_col_types = FALSE
)

if (!"birth_weight_lbs" %in% names(preprocessed_data)) {
  stop("Column birth_weight_lbs is missing from preprocessed_df.csv.")
}

# Plot one overall distribution; participants are not separated or overlaid by
# PCE or any other category.
birthweight_data <- preprocessed_data %>%
  transmute(
    birth_weight_lbs = suppressWarnings(
      as.numeric(as.character(birth_weight_lbs))
    )
  ) %>%
  filter(is.finite(birth_weight_lbs))

birthweight_summary <- birthweight_data %>%
  summarise(
    n = n(),
    mean = mean(birth_weight_lbs),
    sd = stats::sd(birth_weight_lbs),
    median = stats::median(birth_weight_lbs),
    minimum = min(birth_weight_lbs),
    maximum = max(birth_weight_lbs)
  )
print(birthweight_summary)

birthweight_plot <- ggplot(
  birthweight_data,
  aes(x = birth_weight_lbs)
) +
  geom_histogram(
    aes(y = after_stat(density)),
    binwidth = 0.5,
    boundary = 0,
    fill = "#4C8BB6",
    color = "white",
    alpha = 0.8
  ) +
  geom_density(
    color = "#B24C63",
    linewidth = 1
  ) +
  geom_vline(
    xintercept = birthweight_summary$mean,
    color = "#264653",
    linewidth = 0.8,
    linetype = "dashed"
  ) +
  labs(
    title = "Overall Distribution of Birthweight",
    subtitle = paste0(
      "N = ",
      birthweight_summary$n,
      "; dashed line indicates mean = ",
      round(birthweight_summary$mean, 2),
      " lb"
    ),
    x = "Birthweight (lb)",
    y = "Density"
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.title.position = "plot")

output_dir <- file.path(getwd(), "Outputs", "Checks")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

ggsave(
  filename = file.path(output_dir, "birthweight_distribution.png"),
  plot = birthweight_plot,
  width = 9,
  height = 6,
  dpi = 300
)

print(birthweight_plot)

# Unupdated comparison: use only the original pounds component and do not add
# the ounces component.
pregnancy_and_exposures <- read.csv(
  file.path(rawdata_dir, "PregnancyAndExposures.csv"),
  check.names = FALSE
)
original_birthweight_var <- "pex_bm_healthv2_inf_001__02"

if (!original_birthweight_var %in% names(pregnancy_and_exposures)) {
  stop("Original pounds-only birthweight column is missing.")
}

first_numeric_non_missing <- function(x) {
  values <- suppressWarnings(as.numeric(as.character(x)))
  values <- values[is.finite(values)]
  if (length(values) == 0L) NA_real_ else values[[1]]
}

original_birthweight_data <- pregnancy_and_exposures %>%
  arrange(
    participant_id,
    match(
      session_id,
      c("ses-V02", "ses-V01", "ses-V03", "ses-V04", "ses-V05")
    )
  ) %>%
  group_by(participant_id) %>%
  summarise(
    birth_weight_lbs = first_numeric_non_missing(
      .data[[original_birthweight_var]]
    ),
    .groups = "drop"
  ) %>%
  filter(is.finite(birth_weight_lbs))

original_birthweight_summary <- original_birthweight_data %>%
  summarise(
    n = n(),
    mean = mean(birth_weight_lbs),
    sd = stats::sd(birth_weight_lbs),
    median = stats::median(birth_weight_lbs),
    minimum = min(birth_weight_lbs),
    maximum = max(birth_weight_lbs)
  )
print(original_birthweight_summary)

original_pound_counts <- original_birthweight_data %>%
  count(birth_weight_lbs, name = "n_participants") %>%
  arrange(birth_weight_lbs)

original_above_11_summary <- original_birthweight_data %>%
  summarise(
    n_participants_above_11_lbs = sum(birth_weight_lbs > 11),
    percent_above_11_lbs = 100 * mean(birth_weight_lbs > 11)
  )

original_above_11_counts <- original_pound_counts %>%
  filter(birth_weight_lbs > 11)

high_birthweight_gestational_age <- original_birthweight_data %>%
  filter(birth_weight_lbs %in% c(14, 17)) %>%
  left_join(
    preprocessed_data %>%
      select(participant_id, gestational_age),
    by = "participant_id"
  ) %>%
  select(participant_id, birth_weight_lbs, gestational_age) %>%
  arrange(birth_weight_lbs, participant_id)

message("\n========== Counts for each original pounds value ==========")
print(as_tibble(original_pound_counts), n = Inf)
message("\n========== Original pounds values above 11 lb ==========")
print(original_above_11_summary)
print(as_tibble(original_above_11_counts), n = Inf)
message("\n========== Gestational ages for 14 lb and 17 lb records ==========")
print(as_tibble(high_birthweight_gestational_age), n = Inf)

original_birthweight_plot <- ggplot(
  original_birthweight_data,
  aes(x = birth_weight_lbs)
) +
  geom_histogram(
    aes(y = after_stat(density)),
    binwidth = 0.5,
    boundary = 0,
    fill = "#7A9E7E",
    color = "white",
    alpha = 0.8
  ) +
  geom_density(color = "#8C5E58", linewidth = 1) +
  geom_vline(
    xintercept = original_birthweight_summary$mean,
    color = "#264653",
    linewidth = 0.8,
    linetype = "dashed"
  ) +
  labs(
    title = "Original Pounds-Only Birthweight Distribution",
    subtitle = paste0(
      "N = ",
      original_birthweight_summary$n,
      "; dashed line indicates mean = ",
      round(original_birthweight_summary$mean, 2),
      " lb"
    ),
    x = "Original birthweight pounds component",
    y = "Density"
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.title.position = "plot")

ggsave(
  filename = file.path(
    output_dir,
    "birthweight_distribution_original_pounds_only.png"
  ),
  plot = original_birthweight_plot,
  width = 9,
  height = 6,
  dpi = 300
)

print(original_birthweight_plot)
