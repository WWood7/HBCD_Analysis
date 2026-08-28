source("Scripts/00_config.r")
source("Scripts/Helpers/02_eda_helper.r")

library(dplyr)
library(gtsummary)
library(gt)

# Load the merged participant-level dataset.
preprocessed_df <- read.csv(
  file.path(preprocessed_dir, "preprocessed_df.csv"),
  check.names = FALSE
)

# Build and save a separate Table 1 for each treatment definition.
treatment_definitions <- c(
  tox_and_report = "prenatal_cannabis",
  tox_only = "prenatal_cannabis_tox"
)

table1_results <- lapply(names(treatment_definitions), function(definition_name) {
  group_var <- treatment_definitions[[definition_name]]
  table1 <- create_pce_table1(
    preprocessed_df,
    group_var = group_var,
    definition_name = definition_name
  )

  gt::gtsave(
    gtsummary::as_gt(table1),
    filename = file.path(
      eda_output_dir,
      paste0("table1_pce_", definition_name, ".html")
    )
  )

  table1
})
names(table1_results) <- names(treatment_definitions)
