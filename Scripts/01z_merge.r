source("Scripts/00_config.r")
library(dplyr)

# read in the preprocessed tables
demo_ses_table <- read.csv(file.path(preprocessed_dir, "demo_ses_table.csv"))
prenatal_cannabis <- read.csv(file.path(preprocessed_dir, "prenatal_cannabis.csv"))
prenatal_nicotine <- read.csv(file.path(preprocessed_dir, "prenatal_nicotine.csv"))
prenatal_opioid <- read.csv(file.path(preprocessed_dir, "prenatal_opioid.csv"))
prenatal_alcohol <- read.csv(file.path(preprocessed_dir, "prenatal_alcohol.csv"))

# merge the tables
preprocessed_df <- demo_ses_table %>%
    full_join(prenatal_cannabis, by = "participant_id") %>%
    full_join(prenatal_nicotine, by = "participant_id") %>%
    full_join(prenatal_opioid, by = "participant_id") %>%
    full_join(prenatal_alcohol, by = "participant_id")

# write the merged table
write.csv(preprocessed_df, file.path(preprocessed_dir, "preprocessed_df.csv"), row.names = FALSE)