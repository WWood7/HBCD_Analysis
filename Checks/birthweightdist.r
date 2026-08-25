source("Scripts/00_config.r")
source("Scripts/Helpers/01_preprocess_helper.r")

library(dplyr)
library(ggplot2)

# Birth weight is stored in the Pregnancy and Exposures table.
pregnancyandexposures <- read.csv(
    file.path(rawdata_dir, "PregnancyAndExposures.csv"),
    check.names = FALSE
)

birthweight_variable <- "pex_bm_healthv2_inf_001__02"

# Keep one birthweight value per participant, preferring V02 over V01.
birthweight_data <- pregnancyandexposures %>%
  filter(session_id %in% c("ses-V01", "ses-V02")) %>%
  arrange(
    participant_id,
    match(session_id, c("ses-V02", "ses-V01"))
  ) %>%
  group_by(participant_id) %>%
  summarise(
    birth_weight_lbs = first_non_missing(.data[[birthweight_variable]]),
    .groups = "drop"
  ) %>%
  mutate(
    birth_weight_lbs = suppressWarnings(
      as.numeric(as.character(birth_weight_lbs))
    )
  ) %>%
  filter(is.finite(birth_weight_lbs))

birthweight_summary <- birthweight_data %>%
  summarise(
    n = n(),
    mean = mean(birth_weight_lbs),
    sd = sd(birth_weight_lbs),
    median = median(birth_weight_lbs),
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
    title = "Distribution of Birth Weight",
    subtitle = paste0(
      "N = ", birthweight_summary$n,
      "; dashed line indicates mean = ",
      round(birthweight_summary$mean, 2),
      " lb"
    ),
    x = "Birth weight (lb)",
    y = "Density"
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.title.position = "plot")

output_dir <- file.path(getwd(), "Outputs", "Checks")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
birthweight_plot_path <- file.path(
  output_dir,
  "birthweight_distribution.png"
)
ggsave(
  filename = birthweight_plot_path,
  plot = birthweight_plot,
  width = 9,
  height = 6,
  dpi = 300
)
print(birthweight_plot)
