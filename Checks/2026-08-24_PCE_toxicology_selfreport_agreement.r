source("Scripts/00_config.r")
source("Scripts/Helpers/01_preprocess_helper.r")

library(dplyr)
library(ggplot2)

# read in the raw tables
biospecimens <- read.csv(
  file.path(rawdata_dir, "BioSpecimens.csv"),
  check.names = FALSE
)
selfreports <- read.csv(
  file.path(rawdata_dir, "PregnancyAndExposures.csv"),
  check.names = FALSE
)

# merge by both participant_id and session_id
# only keep the records from session V1 and V2
df <- biospecimens %>%
  full_join(selfreports, by = c("participant_id", "session_id")) %>%
  filter(session_id %in% c("ses-V01", "ses-V02"))

# Variables used for toxicology and self-report classifications.
nail_cannabis_variable <-
  "bio_bm_biosample_nails_results_c_delta-9-THC_n_cat"
urine_cannabis_variable <-
  "bio_bm_biosample_urine_results_bio_c_delta-9-THC_u_cat"
assist_variables <- c(
  "pex_bm_assistv1_during__use_003",
  "pex_bm_assistv2_end__use_003"
)

# Return 1 for a positive result, 0 for a negative result, and NA when no
# classifiable result is available.
classify_codes <- function(x, positive_codes, negative_codes) {
  values <- trimws(as.character(x))
  values <- values[!is_missing_value(values)]
  if (any(values %in% as.character(positive_codes))) {
    return(1L)
  }
  if (any(values %in% as.character(negative_codes))) {
    return(0L)
  }
  NA_integer_
}

# Collapse V01 and V02 to one record per participant.
agreement_data <- df %>%
  group_by(participant_id) %>%
  summarise(
    nail_toxicology = classify_codes(
      .data[[nail_cannabis_variable]],
      positive_codes = 1,
      negative_codes = c(0, 4)
    ),
    urine_toxicology = classify_codes(
      .data[[urine_cannabis_variable]],
      positive_codes = 1,
      negative_codes = c(0, 4)
    ),
    assist_self_report = classify_codes(
      unlist(across(all_of(assist_variables)), use.names = FALSE),
      positive_codes = 1,
      negative_codes = 0
    ),
    .groups = "drop"
  ) %>%
  mutate(
    toxicology = case_when(
      nail_toxicology == 1L | urine_toxicology == 1L ~ 1L,
      nail_toxicology == 0L | urine_toxicology == 0L ~ 0L,
      TRUE ~ NA_integer_
    )
  )

safe_ratio <- function(numerator, denominator) {
  if (denominator == 0) NA_real_ else numerator / denominator
}

calculate_agreement <- function(data, toxicology_var, self_report_var, label) {
  paired <- data %>%
    filter(
      !is.na(.data[[toxicology_var]]),
      !is.na(.data[[self_report_var]])
    )

  both_positive <- sum(
    paired[[toxicology_var]] == 1L & paired[[self_report_var]] == 1L
  )
  toxicology_positive_only <- sum(
    paired[[toxicology_var]] == 1L & paired[[self_report_var]] == 0L
  )
  self_report_positive_only <- sum(
    paired[[toxicology_var]] == 0L & paired[[self_report_var]] == 1L
  )
  both_negative <- sum(
    paired[[toxicology_var]] == 0L & paired[[self_report_var]] == 0L
  )
  n_paired <- nrow(paired)
  observed_agreement <- safe_ratio(
    both_positive + both_negative,
    n_paired
  )
  expected_agreement <-
    safe_ratio(
      (both_positive + toxicology_positive_only) *
        (both_positive + self_report_positive_only) +
        (self_report_positive_only + both_negative) *
        (toxicology_positive_only + both_negative),
      n_paired^2
    )

  data.frame(
    comparison = label,
    n_paired = n_paired,
    both_positive = both_positive,
    toxicology_positive_self_report_negative = toxicology_positive_only,
    toxicology_negative_self_report_positive = self_report_positive_only,
    both_negative = both_negative,
    percent_agreement = 100 * observed_agreement,
    positive_agreement = safe_ratio(
      2 * both_positive,
      2 * both_positive + toxicology_positive_only + self_report_positive_only
    ),
    negative_agreement = safe_ratio(
      2 * both_negative,
      2 * both_negative + toxicology_positive_only + self_report_positive_only
    ),
    sensitivity_vs_toxicology = safe_ratio(
      both_positive,
      both_positive + toxicology_positive_only
    ),
    specificity_vs_toxicology = safe_ratio(
      both_negative,
      both_negative + self_report_positive_only
    ),
    cohen_kappa = safe_ratio(
      observed_agreement - expected_agreement,
      1 - expected_agreement
    )
  )
}

comparisons <- data.frame(
  toxicology_var = c(
    "toxicology",
    "nail_toxicology",
    "urine_toxicology"
  ),
  self_report_var = c(
    "assist_self_report",
    "assist_self_report",
    "assist_self_report"
  ),
  label = c(
    "Any toxicology vs ASSIST",
    "Nail toxicology vs ASSIST",
    "Urine toxicology vs ASSIST"
  )
)

agreement_summary <- bind_rows(lapply(seq_len(nrow(comparisons)), function(i) {
  calculate_agreement(
    agreement_data,
    toxicology_var = comparisons$toxicology_var[i],
    self_report_var = comparisons$self_report_var[i],
    label = comparisons$label[i]
  )
}))

classification_summary <- bind_rows(lapply(
  c(
    "toxicology",
    "nail_toxicology",
    "urine_toxicology",
    "assist_self_report"
  ),
  function(variable) {
    agreement_data %>%
      transmute(
        source = variable,
        classification = case_when(
          .data[[variable]] == 1L ~ "Positive",
          .data[[variable]] == 0L ~ "Negative",
          TRUE ~ "Unclassifiable"
        )
      ) %>%
      count(source, classification, name = "n")
  }
))

print(classification_summary)
print(agreement_summary)

# Ribbon plot comparing combined toxicology with ASSIST self-report.
classification_levels <- c("Positive", "Negative")
flow_counts <- agreement_data %>%
  filter(!is.na(toxicology), !is.na(assist_self_report)) %>%
  transmute(
    toxicology_result = if_else(toxicology == 1L, "Positive", "Negative"),
    self_report_result = if_else(
      assist_self_report == 1L,
      "Positive",
      "Negative"
    )
  ) %>%
  count(toxicology_result, self_report_result, name = "n") %>%
  mutate(
    toxicology_result = factor(
      toxicology_result,
      levels = classification_levels
    ),
    self_report_result = factor(
      self_report_result,
      levels = classification_levels
    ),
    flow_id = row_number()
  )

left_strata <- flow_counts %>%
  group_by(toxicology_result) %>%
  summarise(n = sum(n), .groups = "drop") %>%
  arrange(toxicology_result) %>%
  mutate(
    ymin = cumsum(n) - n,
    ymax = cumsum(n),
    axis_x = 0
  )

right_strata <- flow_counts %>%
  group_by(self_report_result) %>%
  summarise(n = sum(n), .groups = "drop") %>%
  arrange(self_report_result) %>%
  mutate(
    ymin = cumsum(n) - n,
    ymax = cumsum(n),
    axis_x = 1
  )

flow_counts <- flow_counts %>%
  left_join(
    left_strata %>%
      select(toxicology_result, left_stratum_ymin = ymin),
    by = "toxicology_result"
  ) %>%
  arrange(toxicology_result, self_report_result) %>%
  group_by(toxicology_result) %>%
  mutate(
    left_ymin = left_stratum_ymin + cumsum(n) - n,
    left_ymax = left_stratum_ymin + cumsum(n)
  ) %>%
  ungroup() %>%
  left_join(
    right_strata %>%
      select(self_report_result, right_stratum_ymin = ymin),
    by = "self_report_result"
  ) %>%
  arrange(self_report_result, toxicology_result) %>%
  group_by(self_report_result) %>%
  mutate(
    right_ymin = right_stratum_ymin + cumsum(n) - n,
    right_ymax = right_stratum_ymin + cumsum(n)
  ) %>%
  ungroup()

# Use smooth-step interpolation to create the alluvial ribbons without an
# additional plotting dependency.
ribbon_data <- bind_rows(lapply(seq_len(nrow(flow_counts)), function(i) {
  flow <- flow_counts[i, ]
  x <- seq(0, 1, length.out = 100)
  smooth_x <- 3 * x^2 - 2 * x^3
  lower <- (1 - smooth_x) * flow$left_ymin + smooth_x * flow$right_ymin
  upper <- (1 - smooth_x) * flow$left_ymax + smooth_x * flow$right_ymax

  data.frame(
    flow_id = flow$flow_id,
    toxicology_result = flow$toxicology_result,
    x = c(x, rev(x)),
    y = c(lower, rev(upper))
  )
}))

strata_data <- bind_rows(
  left_strata %>%
    transmute(
      result = as.character(toxicology_result),
      n,
      ymin,
      ymax,
      axis_x
    ),
  right_strata %>%
    transmute(
      result = as.character(self_report_result),
      n,
      ymin,
      ymax,
      axis_x
    )
) %>%
  mutate(label = paste0(result, "\n(n = ", n, ")"))

agreement_ribbon_plot <- ggplot() +
  geom_polygon(
    data = ribbon_data,
    aes(
      x = x,
      y = y,
      group = flow_id,
      fill = toxicology_result
    ),
    alpha = 0.65,
    color = NA
  ) +
  geom_rect(
    data = strata_data,
    aes(
      xmin = axis_x - 0.035,
      xmax = axis_x + 0.035,
      ymin = ymin,
      ymax = ymax
    ),
    fill = "grey95",
    color = "grey55",
    linewidth = 0.4
  ) +
  geom_text(
    data = strata_data,
    aes(
      x = axis_x,
      y = (ymin + ymax) / 2,
      label = label
    ),
    size = 3.5
  ) +
  scale_x_continuous(
    breaks = c(0, 1),
    labels = c("Combined toxicology", "ASSIST self-report"),
    limits = c(-0.12, 1.12),
    expand = c(0, 0)
  ) +
  scale_fill_manual(
    values = c(
      Positive = "#4C8BB6",
      Negative = "#8FD0C4"
    ),
    drop = FALSE
  ) +
  labs(
    title = "Toxicology Resultsvs Self-Report",
    subtitle = "Participant-level cannabis (delta-9-THC) classifications among participants",
    x = NULL,
    y = "Participants",
    fill = "Toxicology result"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    plot.title.position = "plot"
  )

agreement_plot_dir <- file.path(getwd(), "Outputs", "Checks")
dir.create(agreement_plot_dir, recursive = TRUE, showWarnings = FALSE)
agreement_plot_path <- file.path(
  agreement_plot_dir,
  "pce_toxicology_assist_agreement_ribbon.png"
)
ggsave(
  filename = agreement_plot_path,
  plot = agreement_ribbon_plot,
  width = 10,
  height = 7,
  dpi = 300
)
print(agreement_ribbon_plot)

