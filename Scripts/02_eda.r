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

# Build Table 1 and render it as an HTML table.
table1 <- create_pce_table1(preprocessed_df)
table1_html <- gtsummary::as_gt(table1)

gt::gtsave(
  table1_html,
  filename = file.path(eda_output_dir, "table1_pce.html")
)
